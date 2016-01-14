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
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.
 @ May 17, 2010 - CA
   Code extracted from GReWeightNuXSec and redeveloped in preparation for 
   the Summer 2010 T2K analyses.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNonResonanceBkg::GReWeightNonResonanceBkg() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNonResonanceBkg::~GReWeightNonResonanceBkg()
{

}
//_______________________________________________________________________________________
bool GReWeightNonResonanceBkg::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kXSecTwkDial_RvpCC1pi      ) : 
     case ( kXSecTwkDial_RvpCC2pi      ) : 
     case ( kXSecTwkDial_RvpNC1pi      ) : 
     case ( kXSecTwkDial_RvpNC2pi      ) : 
     case ( kXSecTwkDial_RvnCC1pi      ) : 
     case ( kXSecTwkDial_RvnCC2pi      ) : 
     case ( kXSecTwkDial_RvnNC1pi      ) : 
     case ( kXSecTwkDial_RvnNC2pi      ) : 
     case ( kXSecTwkDial_RvbarpCC1pi   ) : 
     case ( kXSecTwkDial_RvbarpCC2pi   ) :
     case ( kXSecTwkDial_RvbarpNC1pi   ) :
     case ( kXSecTwkDial_RvbarpNC2pi   ) : 
     case ( kXSecTwkDial_RvbarnCC1pi   ) : 
     case ( kXSecTwkDial_RvbarnCC2pi   ) : 
     case ( kXSecTwkDial_RvbarnNC1pi   ) : 
     case ( kXSecTwkDial_RvbarnNC2pi   ) : 

          handle = true;
          break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightNonResonanceBkg::SetSystematic(GSyst_t syst, double twk_dial)
{
   switch(syst) {
     case ( kXSecTwkDial_RvpCC1pi      ) : 
     case ( kXSecTwkDial_RvpCC2pi      ) : 
     case ( kXSecTwkDial_RvpNC1pi      ) : 
     case ( kXSecTwkDial_RvpNC2pi      ) : 
     case ( kXSecTwkDial_RvnCC1pi      ) : 
     case ( kXSecTwkDial_RvnCC2pi      ) : 
     case ( kXSecTwkDial_RvnNC1pi      ) : 
     case ( kXSecTwkDial_RvnNC2pi      ) : 
     case ( kXSecTwkDial_RvbarpCC1pi   ) : 
     case ( kXSecTwkDial_RvbarpCC2pi   ) :
     case ( kXSecTwkDial_RvbarpNC1pi   ) :
     case ( kXSecTwkDial_RvbarpNC2pi   ) : 
     case ( kXSecTwkDial_RvbarnCC1pi   ) : 
     case ( kXSecTwkDial_RvbarnCC2pi   ) : 
     case ( kXSecTwkDial_RvbarnNC1pi   ) : 
     case ( kXSecTwkDial_RvbarnNC2pi   ) : 

          fRTwkDial[syst] = twk_dial;
          break;

     default:
          break;
   }
}
//_______________________________________________________________________________________
void GReWeightNonResonanceBkg::Reset(void)
{
  map<GSyst_t,double>::iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    it->second = 0.;
  }

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNonResonanceBkg::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  map<GSyst_t,double>::iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    GSyst_t syst = it->first;
    double twk_dial = it->second;
    double curr = fRDef[syst] * (1 + twk_dial * fracerr->OneSigmaErr(syst));
    fRCurr[syst] = TMath::Max(0.,curr);
  }
}
//_______________________________________________________________________________________
double GReWeightNonResonanceBkg::CalcWeight(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  bool is_dis = interaction->ProcInfo().IsDeepInelastic();
  if(!is_dis) return 1.;

  bool selected = true;
  double W = interaction->Kine().W(selected);
  bool in_transition = (W<fWmin);
  if(!in_transition) return 1.;

  int probe  = interaction->InitState().ProbePdg();
  int hitnuc = interaction->InitState().Tgt().HitNucPdg();
  InteractionType_t itype = interaction->ProcInfo().InteractionTypeId();

  int nhadmult = 0;
  int nnuc     = 0;
  int npi      = 0;

  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     int pdgc = p->Pdg();        
     int imom = p->FirstMother();
     if(imom == -1) continue;

     GHepParticle * mom = event.Particle(imom);
     if(!mom) continue;

     if(mom->Pdg() == kPdgHadronicSyst) 
     {
        nhadmult++;
        if ( pdg::IsNucleon(pdgc) ) { nnuc++; }
        if ( pdg::IsPion   (pdgc) ) { npi++;  }
     }
  }//p

  if(nhadmult < 2 || nhadmult > 3) return 1.;
  if(nnuc != 1) return 1.;

  GSyst_t syst = GSyst::RBkg(itype, probe, hitnuc, npi);
  if(syst == kNullSystematic) return 1.;

  double curr = fRCurr[syst];
  double def  = fRDef[syst];

  if(def>0. && curr>=0.) {
    double wght = curr / def;
    return wght;
  }

  return 1.;
}
//_______________________________________________________________________________________
double GReWeightNonResonanceBkg::CalcChisq()
{
  double chisq = 0.;
  map<GSyst_t,double>::const_iterator it = fRTwkDial.begin();
  for( ; it != fRTwkDial.end(); ++it) {
    double twk_dial = it->second;
    chisq += TMath::Power(twk_dial, 2.);
  }
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNonResonanceBkg::Init(void)
{
  this->SetWminCut(2.0*units::GeV);

  // Get the default  parameters 
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * user_config = conf_pool->GlobalParameterList();

  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvpCC1pi,    user_config->GetDouble("DIS-HMultWgt-vp-CC-m2" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvpCC2pi,    user_config->GetDouble("DIS-HMultWgt-vp-CC-m3" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvpNC1pi,    user_config->GetDouble("DIS-HMultWgt-vp-NC-m2" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvpNC2pi,    user_config->GetDouble("DIS-HMultWgt-vp-NC-m3" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvnCC1pi,    user_config->GetDouble("DIS-HMultWgt-vn-CC-m2" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvnCC2pi,    user_config->GetDouble("DIS-HMultWgt-vn-CC-m3" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvnNC1pi,    user_config->GetDouble("DIS-HMultWgt-vn-NC-m2" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvnNC2pi,    user_config->GetDouble("DIS-HMultWgt-vn-NC-m3" )));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarpCC1pi, user_config->GetDouble("DIS-HMultWgt-vbp-CC-m2")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarpCC2pi, user_config->GetDouble("DIS-HMultWgt-vbp-CC-m3")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarpNC1pi, user_config->GetDouble("DIS-HMultWgt-vbp-NC-m2")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarpNC2pi, user_config->GetDouble("DIS-HMultWgt-vbp-NC-m3")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarnCC1pi, user_config->GetDouble("DIS-HMultWgt-vbn-CC-m2")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarnCC2pi, user_config->GetDouble("DIS-HMultWgt-vbn-CC-m3")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarnNC1pi, user_config->GetDouble("DIS-HMultWgt-vbn-NC-m2")));
  fRDef.insert(map<GSyst_t,double>::value_type(
     kXSecTwkDial_RvbarnNC2pi, user_config->GetDouble("DIS-HMultWgt-vbn-NC-m3")));

  map<GSyst_t,double>::const_iterator it = fRDef.begin();
  for( ; it != fRDef.end(); ++it) {
    fRCurr.insert(map<GSyst_t,double>::value_type(it->first, it->second));
  }
}
//_______________________________________________________________________________________
