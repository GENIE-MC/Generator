//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
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
 @ Oct 20, 2010 - CA
   Make static consts kModeABCV12u and kModeABCV12uShape public so as to
   aid external configuration.
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>
#include <TFile.h>
#include <TNtupleD.h>

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

//#define _G_REWEIGHT_DIS_DEBUG_

using namespace genie;
using namespace genie::rew;
   
const int GReWeightNuXSecDIS::kModeABCV12u;
const int GReWeightNuXSecDIS::kModeABCV12uShape;

//_______________________________________________________________________________________
GReWeightNuXSecDIS::GReWeightNuXSecDIS() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecDIS::~GReWeightNuXSecDIS()
{
#ifdef _G_REWEIGHT_DIS_DEBUG_   
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecDIS::IsHandled(GSyst_t syst)
{
   bool handle = false;

   switch(syst) {
     case ( kXSecTwkDial_AhtBYshape  ) : 
     case ( kXSecTwkDial_BhtBYshape  ) : 
     case ( kXSecTwkDial_CV1uBYshape ) : 
     case ( kXSecTwkDial_CV2uBYshape ) : 

       if      (fMode == kModeABCV12u     ) handle = false;
       else if (fMode == kModeABCV12uShape) handle = true;
       break;

     case ( kXSecTwkDial_AhtBY  ) : 
     case ( kXSecTwkDial_BhtBY  ) : 
     case ( kXSecTwkDial_CV1uBY ) : 
     case ( kXSecTwkDial_CV2uBY ) : 

       if      (fMode == kModeABCV12u     ) handle = true;
       else if (fMode == kModeABCV12uShape) handle = false;
       break;

     default:
       handle = false;
       break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSecDIS::SetSystematic(GSyst_t syst, double twk_dial)
{
   if(! this->IsHandled(syst)) return;

   switch(syst) {
     case ( kXSecTwkDial_AhtBY       ) : 
     case ( kXSecTwkDial_AhtBYshape  ) : 
       fAhtBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_BhtBY       ) : 
     case ( kXSecTwkDial_BhtBYshape  ) : 
       fBhtBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_CV1uBY      ) : 
     case ( kXSecTwkDial_CV1uBYshape ) : 
       fCV1uBYTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_CV2uBY      ) : 
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

  double fracerr_aht  = 0.;
  double fracerr_bht  = 0.;
  double fracerr_cv1u = 0.;
  double fracerr_cv2u = 0.;

  if(fMode == kModeABCV12u) {
      fracerr_aht  = fracerr->OneSigmaErr (kXSecTwkDial_AhtBY );
      fracerr_bht  = fracerr->OneSigmaErr (kXSecTwkDial_BhtBY );
      fracerr_cv1u = fracerr->OneSigmaErr (kXSecTwkDial_CV1uBY);
      fracerr_cv2u = fracerr->OneSigmaErr (kXSecTwkDial_CV2uBY);
  }
  else
  if(fMode == kModeABCV12uShape) {
      fracerr_aht  = fracerr->OneSigmaErr (kXSecTwkDial_AhtBYshape );
      fracerr_bht  = fracerr->OneSigmaErr (kXSecTwkDial_BhtBYshape );
      fracerr_cv1u = fracerr->OneSigmaErr (kXSecTwkDial_CV1uBYshape);
      fracerr_cv2u = fracerr->OneSigmaErr (kXSecTwkDial_CV2uBYshape);
  }

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
  bool tweaked = 
      (TMath::Abs(fAhtBYTwkDial)  > controls::kASmallNum) ||
      (TMath::Abs(fBhtBYTwkDial)  > controls::kASmallNum) ||
      (TMath::Abs(fCV1uBYTwkDial) > controls::kASmallNum) ||
      (TMath::Abs(fCV2uBYTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  bool is_dis = interaction->ProcInfo().IsDeepInelastic();
  if(!is_dis) return 1.;

  bool charm = event.Summary()->ExclTag().IsCharmEvent(); // skip DIS charm
  if(charm) return 1.;

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

  double wght = 1.;

  if(fMode == kModeABCV12u) {
     wght = this->CalcWeightABCV12u(event);
  }
  else
  if(fMode == kModeABCV12uShape) {
     wght = this->CalcWeightABCV12uShape(event);
  }

#ifdef _G_REWEIGHT_DIS_DEBUG_
  double E = interaction->InitState().ProbeE(kRfHitNucRest);
  double x = interaction->Kine().x(true);
  double y = interaction->Kine().y(true);
  int ccnc = (is_cc) ? 1 : 0;
  int nuc  = interaction->InitState().Tgt().HitNucPdg();
  int qrk  = interaction->InitState().Tgt().HitQrkPdg();
  int sea  = (interaction->InitState().Tgt().HitSeaQrk()) ? 1 : 0;
  fTestNtp->Fill(E,x,y,nupdg,nuc,qrk,sea,ccnc,wght);
#endif

  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSecDIS::CalcWeightABCV12u(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double twk_xsec   = fXSecModel->XSec(interaction, kPSxyfE);
  double weight = old_weight * (twk_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecDIS::CalcWeightABCV12uShape(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double twk_xsec   = fXSecModel->XSec(interaction, kPSxyfE);
  double weight = old_weight * (twk_xsec/old_xsec);

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double twk_integrated_xsec = fXSecModel    -> Integral(interaction);   

  assert(twk_integrated_xsec > 0);
  weight *= (old_integrated_xsec/twk_integrated_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecDIS::CalcChisq(void)
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

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->SetMode   (kModeABCV12u);

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

#ifdef _G_REWEIGHT_DIS_DEBUG_
  fTestFile = new TFile("./dis_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:x:y:nu:nuc:qrk:sea:ccnc:wght");
#endif
}
//_______________________________________________________________________________________


