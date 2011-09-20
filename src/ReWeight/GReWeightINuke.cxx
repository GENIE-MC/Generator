//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Jul 29, 2011 - SD,AM
   Mean free path is now function of Z too.

*/
//____________________________________________________________________________
#define _T2KRW_REWEIGHT_INUKE_DEBUG_NTP_
#define _G_REWEIGHT_INUKE_DEBUG_NTP_
#include <cassert>

#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightUtils.h"
#include "ReWeight/GSystUncertainty.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightINuke::GReWeightINuke() :
GReWeightI()
{
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
  fTestNtp = new TNtuple("testntp","","pdg:E:tweakDial:d:dmfp:fate:interact");
#endif
}
//_______________________________________________________________________________________
GReWeightINuke::~GReWeightINuke()
{
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
  TFile f("./inuke_reweight_test.root","recreate");
  fTestNtp->Write();
  f.Close();
  delete fTestNtp;
#endif
}
//_______________________________________________________________________________________
bool GReWeightINuke::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kINukeTwkDial_MFP_pi      ) :
     case ( kINukeTwkDial_MFP_N       ) :
     case ( kINukeTwkDial_FrCEx_pi    ) :
     case ( kINukeTwkDial_FrElas_pi   ) :
     case ( kINukeTwkDial_FrInel_pi   ) :
     case ( kINukeTwkDial_FrAbs_pi    ) :
     case ( kINukeTwkDial_FrPiProd_pi ) :
     case ( kINukeTwkDial_FrCEx_N     ) :
     case ( kINukeTwkDial_FrElas_N    ) :
     case ( kINukeTwkDial_FrInel_N    ) :
     case ( kINukeTwkDial_FrAbs_N     ) :
     case ( kINukeTwkDial_FrPiProd_N  ) :
          handle = true;
          break;

     default:
          handle = false;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightINuke::SetSystematic(GSyst_t syst, double val)
{
  if(this->IsHandled(syst)) {
     fINukeRwParams.SetCurTwkDial(syst, val);
  }
}
//_______________________________________________________________________________________
void GReWeightINuke::Reset(void)
{
  fINukeRwParams.Reset();
  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightINuke::Reconfigure(void)
{
  fINukeRwParams.Reconfigure();  
}
//_______________________________________________________________________________________
double GReWeightINuke::CalcWeight(const EventRecord & event) 
{  
  // get the atomic mass number for the hit nucleus
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.0;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.0;
  if (Z<=1) return 1.0;

  double event_weight  = 1.0;

  // Loop over stdhep entries and only calculate weights for particles. 
  // All particles that are not hadrons generated inside the nucleus are given weights of 1.0
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     // Skip particles not rescattered by the actual hadron transport code      
     int  pdgc       = p->Pdg();
     bool is_pion    = pdg::IsPion   (pdgc);
     bool is_nucleon = pdg::IsNucleon(pdgc);      
     if(!is_pion && !is_nucleon) 
     { 
        continue; 
     }

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus) 
     {
        continue; 
     }

     // Get 4-momentum and 4-position
     TLorentzVector x4 (p->Vx(), p->Vy(), p->Vz(), 0.    );
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());

     // Determine the interaction type for current hadron in nucleus, if any
     bool interacted = (p->RescatterCode() != -1);
     INukeFateHA_t hadron_fate = (interacted) ?
          (INukeFateHA_t) p->RescatterCode() : kIHAFtUndefined;
    
     LOG("ReW", pDEBUG) 
        << "Reweighting hadron: PDG code = " << pdgc 
        << ", Fate id = "  << hadron_fate 
        << " (" << INukeHadroFates::AsString(hadron_fate) << ")";

     // Init current hadron weights
     double w_mfp  = 1.0;
     double w_fate = 1.0;

     // Check which weights need to be calculated (only if relevant params were tweaked)
     bool calc_w_mfp  = fINukeRwParams.MeanFreePathParams(pdgc)->IsTweaked();
     bool calc_w_fate = fINukeRwParams.FateParams(pdgc)->IsTweaked();
       
     // Compute weight to account for changes in the total rescattering probability
     if(calc_w_mfp)
     {
        double mfp_scale_factor = fINukeRwParams.MeanFreePathParams(pdgc)->ScaleFactor();
        w_mfp = utils::rew::MeanFreePathWeight(pdgc,x4,p4,A,Z,mfp_scale_factor,interacted);

        // Debug info
#ifdef _T2KRW_REWEIGHT_INUKE_DEBUG_NTP_
        double d        = utils::intranuke::Dist2Exit(x4,p4,A);
        double d_mfp    = utils::intranuke::Dist2ExitMFP(pdgc,x4,p4,A);
        double Eh       = p->E();
	double iflag    = (interacted) ? 1 : -1;
        fTestNtp->Fill(pdgc, Eh, mfp_scale_factor, d, d_mfp, hadron_fate, iflag);
#endif
     } // calculate mfp weight?

     // Compute weight to account for changes in relative fractions of reaction channels
     if(calc_w_fate && interacted)
     {    
        // JIMTODO - Need to deal with normalisation properly
        double fate_fraction_scale_factor = 
             fINukeRwParams.FateParams(pdgc)->ScaleFactor(
                  GSyst::INukeFate2GSyst(hadron_fate,pdgc), p4);
        w_fate = fate_fraction_scale_factor;
     } 

     // Calculate the current hadron weight
     double hadron_weight = w_mfp * w_fate;

     LOG("ReW", pNOTICE) 
        << "Reweighted hadron: PDG code = " << pdgc 
        << ", Fate id = "  << hadron_fate 
        << " (" << INukeHadroFates::AsString(hadron_fate) << ") >"
        << " w_mfp = "  << w_mfp
        <<", w_fate = " << w_fate;

     // Update the current event weight
     event_weight *= hadron_weight;
     
  }//particle loop
  
  return event_weight;
}
//_______________________________________________________________________________________
double GReWeightINuke::CalcChisq(void)
{
  return fINukeRwParams.ChisqPenalty();
}
//_______________________________________________________________________________________
