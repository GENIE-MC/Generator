//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Jul 29, 2011 - SD,AM
   Mean free path is now function of Z too.
 @ Sep 20, 2011 - CA
   Determine 'no interaction' by looking-up the FSI code set by INTRANUKE.
   INTRANUKE now sets a value in all cases rather than leaving it unset for
   particles which escape the target nucleus.

*/
//____________________________________________________________________________

//#define _G_REWEIGHT_INUKE_DEBUG_NTP_

#include <cassert>
#include <cstdlib>

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
  fTestFile = new TFile("./intranuke_reweight_test.root","recreate");
  fTestNtp  = new TNtuple("testntp","","pdg:E:mfp_twk_dial:d:d_mfp:fate:interact:w_mfp:w_fate");
#endif
}
//_______________________________________________________________________________________
GReWeightINuke::~GReWeightINuke()
{
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
  assert(fTestFile);
  assert(fTestNtp);
  fTestFile->cd();
  fTestNtp->Write();
  fTestFile->Close();
  delete fTestFile;
  //delete fTestNtp;
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
     fINukeRwParams.SetTwkDial(syst, val);
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
  int ip=-1;
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
     ip++;

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

     // Determine the interaction type for current hadron in nucleus, if any
     int fsi_code = p->RescatterCode();
     LOG("ReW", pDEBUG) 
        << "Attempting to reweight hadron at position = " << ip 
        << " with PDG code = " << pdgc 
        << " and FSI code = "  << fsi_code 
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ")";
     if(fsi_code == -1 || fsi_code == (int)kIHAFtUndefined) {
       LOG("ReW", pFATAL) << "INTRANUKE didn't set a valid rescattering code for event in position: " << ip;
       LOG("ReW", pFATAL) << "Here is the problematic event:";
       LOG("ReW", pFATAL) << event;
       exit(1);
     }
     bool escaped    = (fsi_code == (int)kIHAFtNoInteraction);
     bool interacted = !escaped;

     // Get 4-momentum and 4-position
     TLorentzVector x4 (p->Vx(), p->Vy(), p->Vz(), 0.    );
     TLorentzVector p4 (p->Px(), p->Py(), p->Pz(), p->E());

     // Init current hadron weights
     double w_mfp  = 1.0;
     double w_fate = 1.0;

     // Check which weights need to be calculated (only if relevant params were tweaked)
     bool calc_w_mfp  = fINukeRwParams.MeanFreePathParams(pdgc)->IsTweaked();
     bool calc_w_fate = fINukeRwParams.FateParams(pdgc)->IsTweaked();

     // Compute weight to account for changes in the total rescattering probability
     double mfp_scale_factor = 1.;
     if(calc_w_mfp)
     {
        mfp_scale_factor = fINukeRwParams.MeanFreePathParams(pdgc)->ScaleFactor();
        w_mfp = utils::rew::MeanFreePathWeight(pdgc,x4,p4,A,Z,mfp_scale_factor,interacted);
     } // calculate mfp weight?

     // Compute weight to account for changes in relative fractions of reaction channels
     if(calc_w_fate && interacted)
     {    
        double fate_fraction_scale_factor = 
             fINukeRwParams.FateParams(pdgc)->ScaleFactor(
                  GSyst::INukeFate2GSyst((INukeFateHA_t)fsi_code,pdgc), p4);
        w_fate = fate_fraction_scale_factor;
     } 

     // Calculate the current hadron weight
     double hadron_weight = w_mfp * w_fate;

     LOG("ReW", pNOTICE) 
        << "Reweighted hadron at position = " << ip
        << " with PDG code = " << pdgc 
        << ", FSI code = "  << fsi_code 
        << " (" << INukeHadroFates::AsString((INukeFateHA_t)fsi_code) << ") :"
        << " w_mfp = "  << w_mfp
        <<", w_fate = " << w_fate;

     // Debug info
#ifdef _G_REWEIGHT_INUKE_DEBUG_NTP_
     double d        = utils::intranuke::Dist2Exit(x4,p4,A);
     double d_mfp    = utils::intranuke::Dist2ExitMFP(pdgc,x4,p4,A,Z);
     double Eh       = p->E();
     double iflag    = (interacted) ? 1 : 0;
     fTestNtp->Fill(pdgc, Eh, mfp_scale_factor, d, d_mfp, fsi_code, iflag, w_mfp, w_fate);
#endif

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
