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
 @ Oct 20, 2010 - CA
   Made static consts `kModeMaMv' and `kModeNormAndMaMvShape' public to
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
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"
#include "Registry/Registry.h"

//#define _G_REWEIGHT_CCRES_DEBUG_

using namespace genie;
using namespace genie::rew;

const int GReWeightNuXSecCCRES::kModeMaMv;
const int GReWeightNuXSecCCRES::kModeNormAndMaMvShape;

//_______________________________________________________________________________________
GReWeightNuXSecCCRES::GReWeightNuXSecCCRES() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCRES::~GReWeightNuXSecCCRES()
{
#ifdef _G_REWEIGHT_CCRES_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCRES::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {

     case ( kXSecTwkDial_NormCCRES    ) :
     case ( kXSecTwkDial_MaCCRESshape ) :
     case ( kXSecTwkDial_MvCCRESshape ) :
       if(fMode==kModeNormAndMaMvShape) { 
          handle = true;  
       } else { 
          handle = false; 
       }
       break;

     case ( kXSecTwkDial_MaCCRES ) :
     case ( kXSecTwkDial_MvCCRES ) :
       if(fMode==kModeMaMv) { 
          handle = true;  
       } else { 
          handle = false; 
       }
       break;

     default:
          handle = false;
          break;
   }

   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCRES::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_NormCCRES ) :
      fNormTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_MaCCRESshape ) :
    case ( kXSecTwkDial_MaCCRES ) :
      fMaTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_MvCCRESshape ) :
    case ( kXSecTwkDial_MvCCRES ) :
      fMvTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCRES::Reset(void)
{
  fNormTwkDial = 0.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvCurr      = fMvDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCRES::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  if(fMode==kModeMaMv) {   
     double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCCRES);
     double fracerr_mv = fracerr->OneSigmaErr(kXSecTwkDial_MvCCRES);
     fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
     fMvCurr = fMvDef * (1. + fMvTwkDial * fracerr_mv);
  }
  else
  if(fMode==kModeNormAndMaMvShape) { 
     double fracerr_norm = fracerr->OneSigmaErr(kXSecTwkDial_NormCCRES);
     double fracerr_mash = fracerr->OneSigmaErr(kXSecTwkDial_MaCCRESshape);
     double fracerr_mvsh = fracerr->OneSigmaErr(kXSecTwkDial_MvCCRESshape);
     fNormCurr = fNormDef * (1. + fNormTwkDial * fracerr_norm);
     fMaCurr   = fMaDef   * (1. + fMaTwkDial   * fracerr_mash);
     fMvCurr   = fMvDef   * (1. + fMvTwkDial   * fracerr_mvsh);
  }

  fNormCurr = TMath::Max(0., fNormCurr);
  fMaCurr   = TMath::Max(0., fMaCurr  );
  fMvCurr   = TMath::Max(0., fMvCurr  );

  Registry & r = *fXSecModelConfig;

  r.Set(fMaPath, fMaCurr); 
  r.Set(fMvPath, fMvCurr); 
  fXSecModel->Configure(r);

//LOG("ReW, pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCRES::CalcWeight(const genie::EventRecord & event) 
{
  bool is_res = event.Summary()->ProcInfo().IsResonant();
  bool is_cc  = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_res || !is_cc) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  if(fMode==kModeMaMv) {
     double wght = this->CalcWeightMaMv(event);
     return wght;
  } 
  else 
  if(fMode==kModeNormAndMaMvShape) {
     double wght = 
         this->CalcWeightNorm      (event) *
         this->CalcWeightMaMvShape (event);
     return wght;
  }

  return 1.;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCRES::CalcChisq()
{
  double chisq = 0.;
  if(fMode==kModeMaMv) {   
     chisq += TMath::Power(fMaTwkDial, 2.);
     chisq += TMath::Power(fMvTwkDial, 2.);
  }
  else
  if(fMode==kModeNormAndMaMvShape) { 
     chisq += TMath::Power(fNormTwkDial, 2.);
     chisq += TMath::Power(fMaTwkDial,   2.);
     chisq += TMath::Power(fMvTwkDial,   2.);
  }
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCRES::Init(void)
{
  AlgId id("genie::ReinSehgalRESPXSec","Default");

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->SetMode(kModeNormAndMaMvShape);

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  this->SetMaPath("Ma");
  this->SetMvPath("Mv");

  fNormTwkDial = 0.;
  fNormDef     = 1.;
  fNormCurr    = fNormDef;
  fMaTwkDial   = 0.; 
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fMvTwkDial   = 0.; 
  fMvDef       = fXSecModelConfig->GetDouble(fMvPath);
  fMvCurr      = fMvDef;

#ifdef _G_REWEIGHT_CCRES_DEBUG_
  fTestFile = new TFile("./ccres_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:W:wght");
#endif
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCRES::CalcWeightNorm(const genie::EventRecord & /*event*/) 
{
  bool tweaked = (TMath::Abs(fNormTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fNormCurr;
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCRES::CalcWeightMaMv(const genie::EventRecord & event) 
{
  bool tweaked = 
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSWQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCRES::CalcWeightMaMvShape(const genie::EventRecord & event) 
{
  bool tweaked = 
     (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
     (TMath::Abs(fMvTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction * interaction = event.Summary();

  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSWQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

//double old_integrated_xsec = event.XSec();
  double old_integrated_xsec = fXSecModelDef -> Integral(interaction);
  double twk_integrated_xsec = fXSecModel    -> Integral(interaction);
  assert(twk_integrated_xsec > 0);
  new_weight *= (old_integrated_xsec/twk_integrated_xsec);

//LOG("ReW", pDEBUG) << "integrated cross section (old) = " << old_integrated_xsec;
//LOG("ReW", pDEBUG) << "integrated cross section (twk) = " << twk_integrated_xsec;
//LOG("ReW", pDEBUG) << "new weight (normalized to const integral) = " << new_weight;

  interaction->KinePtr()->ClearRunningValues();

#ifdef _G_REWEIGHT_CCRES_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  double W  = interaction->Kine().W(true);
  fTestNtp->Fill(E,Q2,W,new_weight);
#endif

  return new_weight;
}
//_______________________________________________________________________________________
