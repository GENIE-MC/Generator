//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

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
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"
#include "Registry/Registry.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecDIS::GReWeightNuXSecDIS() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecDIS::~GReWeightNuXSecDIS()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecDIS::IsHandled(GSyst_t syst)
{
   switch(syst) {
     case ( kXSecTwkDial_AhtBYshape  ) : 
     case ( kXSecTwkDial_BhtBYshape  ) : 
     case ( kXSecTwkDial_CV1uBYshape ) : 
     case ( kXSecTwkDial_CV2uBYshape ) : 
       return true;
       break;
     default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecDIS::SetSystematic(GSyst_t syst, double twk_dial)
{
   switch(syst) {
     case ( kXSecTwkDial_AhtBYshape  ) : 
       fAhtBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_BhtBYshape  ) : 
       fBhtBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_CV1uBYshape ) : 
       fCV1uBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_CV2uBYshape ) : 
       fCV2uBYTwkDial = twk_dial;
       break;
     default:
       return;
       break;
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSecDIS::Reset(void)
{
  fAhtBYTwkDial  = 0.;        
  fBhtBYTwkDial  = 0.;        
  fCV1uBYTwkDial = 0.;      
  fCV2uBYTwkDial = 0.;       

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecDIS::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  double fracerr_aht  = fracerr->OneSigmaErr (kXSecTwkDial_AhtBYshape );
  double fracerr_bht  = fracerr->OneSigmaErr (kXSecTwkDial_BhtBYshape );
  double fracerr_cv1u = fracerr->OneSigmaErr (kXSecTwkDial_CV1uBYshape);
  double fracerr_cv2u = fracerr->OneSigmaErr (kXSecTwkDial_CV2uBYshape);

  fAhtBYCur  = fAhtBYDef  * (1. + fAhtBYTwkDial  * fracerr_aht );
  fBhtBYCur  = fBhtBYDef  * (1. + fBhtBYTwkDial  * fracerr_bht );
  fCV1uBYCur = fCV1uBYDef * (1. + fCV1uBYTwkDial * fracerr_cv1u);
  fCV2uBYCur = fCV2uBYDef * (1. + fCV2uBYTwkDial * fracerr_cv2u);

  Registry & r = *fXSecModelConfig;

  r.Set(fAhtBYPath,  fAhtBYCur ); 
  r.Set(fBhtBYPath,  fBhtBYCur ); 
  r.Set(fCV1uBYPath, fCV1uBYCur); 
  r.Set(fCV2uBYPath, fCV2uBYCur); 

  fXSecModel->Configure(r);

//LOG("ReW", pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecDIS::CalcWeight(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  bool is_dis = interaction->ProcInfo().IsDeepInelastic();
  if(!is_dis) return 1.;

  bool charm = event.Summary()->ExclTag().IsCharmEvent(); // skip DIS charm
  if(charm) return 1.;

  bool tweaked = 
      (TMath::Abs(fAhtBYTwkDial)  > controls::kASmallNum) ||
      (TMath::Abs(fBhtBYTwkDial)  > controls::kASmallNum) ||
      (TMath::Abs(fCV1uBYTwkDial) > controls::kASmallNum) ||
      (TMath::Abs(fCV2uBYTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  bool is_cc  = interaction->ProcInfo().IsWeakCC();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  bool selected = true;
  double W  = interaction->Kine().W (selected);
  double Q2 = interaction->Kine().Q2(selected);
  bool passes_kine_cuts = (W>=fWmin && Q2>=fQ2min);
  if(!passes_kine_cuts) return 1.;

  //
  // calculate weight
  //

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSxyfE);
  double new_weight = old_weight * (new_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecDIS::CalcChisq()
{
  double chisq = 
      TMath::Power(fAhtBYTwkDial,  2.) +
      TMath::Power(fBhtBYTwkDial,  2.) +
      TMath::Power(fCV1uBYTwkDial, 2.) +
      TMath::Power(fCV2uBYTwkDial, 2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecDIS::Init(void)
{
  AlgId id("genie::QPMDISPXSec","Default");

  AlgFactory * algf = AlgFactory::Instance();
  Algorithm * alg = algf->AdoptAlgorithm(id);

  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  this->SetWminCut (1.7*units::GeV);
  this->SetQ2minCut(1.0*units::GeV*units::GeV);

  this->SetAhtBYPath ("SFAlg/A");
  this->SetBhtBYPath ("SFAlg/B");
  this->SetCV1uBYPath("SFAlg/Cv1U");
  this->SetCV2uBYPath("SFAlg/Cv2U");

  fAhtBYTwkDial  = 0;   
  fBhtBYTwkDial  = 0;   
  fCV1uBYTwkDial = 0;  
  fCV2uBYTwkDial = 0;  
  fAhtBYDef  = fXSecModelConfig->GetDouble(fAhtBYPath);       
  fBhtBYDef  = fXSecModelConfig->GetDouble(fBhtBYPath);
  fCV1uBYDef = fXSecModelConfig->GetDouble(fCV1uBYPath);
  fCV2uBYDef = fXSecModelConfig->GetDouble(fCV2uBYPath);
  fAhtBYCur  = fAhtBYDef;       
  fBhtBYCur  = fBhtBYDef;       
  fCV1uBYCur = fCV1uBYDef;      
  fCV2uBYCur = fCV2uBYDef;      
}
//_______________________________________________________________________________________


