//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 25, 2010 - CA
   First included in v2.7.1. Added option to reweight NCEL axial mass M_{A} 
   and the strange axial form factor eta.
*/
//____________________________________________________________________________

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
#include "ReWeight/GReWeightNuXSecNCEL.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"
#include "ReWeight/GReWeightUtils.h"
#include "Registry/Registry.h"

//#define _G_REWEIGHT_NCEL_DEBUG_

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecNCEL::GReWeightNuXSecNCEL() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecNCEL::~GReWeightNuXSecNCEL()
{
#ifdef _G_REWEIGHT_NCEL_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecNCEL::IsHandled(GSyst_t syst)
{
   bool handle;
   switch(syst) {
     case ( kXSecTwkDial_MaNCEL  ) : handle = true; break;
     case ( kXSecTwkDial_EtaNCEL ) : handle = true; break;
     default:
          handle = false;
          break;
   }
   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_MaNCEL ) :
      fMaTwkDial = twk_dial;
      break;
    case ( kXSecTwkDial_EtaNCEL ) :
      fEtaTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Reset(void)
{
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;
  fEtaTwkDial  = 0.; 
  fEtaCurr     = fEtaDef;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  int sign_matwk  = utils::rew::Sign(fMaTwkDial );
  int sign_etatwk = utils::rew::Sign(fEtaTwkDial);
  
  double fracerr_ma  = fracerr->OneSigmaErr(kXSecTwkDial_MaNCEL,  sign_matwk );
  double fracerr_eta = fracerr->OneSigmaErr(kXSecTwkDial_EtaNCEL, sign_etatwk);
  
  fMaCurr  = fMaDef  * (1. + fMaTwkDial  * fracerr_ma);
  fEtaCurr = fEtaDef * (1. + fEtaTwkDial * fracerr_eta);

  fMaCurr  = TMath::Max(0., fMaCurr  );
  fEtaCurr = TMath::Max(0., fEtaCurr );

  Registry & r = *fXSecModelConfig;

  r.Set(fMaPath,  fMaCurr ); 
  r.Set(fEtaPath, fEtaCurr); 
  fXSecModel->Configure(r);

//LOG("ReW, pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecNCEL::CalcWeight(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  bool is_qe = interaction->ProcInfo().IsQuasiElastic();
  bool is_nc = interaction->ProcInfo().IsWeakNC();
  if(!is_qe || !is_nc) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  bool tweaked_ma  = (TMath::Abs(fMaTwkDial ) > controls::kASmallNum);
  bool tweaked_eta = (TMath::Abs(fEtaTwkDial) > controls::kASmallNum);
  bool tweaked     = tweaked_ma || tweaked_eta;
  if(!tweaked) return 1.0;

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double old_xsec   = event.DiffXSec();
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, kPSQ2fE);
  double new_weight = old_weight * (new_xsec/old_xsec);

//LOG("ReW", pDEBUG) << "differential cross section (old) = " << old_xsec;
//LOG("ReW", pDEBUG) << "differential cross section (new) = " << new_xsec;
//LOG("ReW", pDEBUG) << "event generation weight = " << old_weight;
//LOG("ReW", pDEBUG) << "new weight = " << new_weight;

#ifdef _G_REWEIGHT_NCEL_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(E,Q2,new_weight);
#endif

  interaction->KinePtr()->ClearRunningValues();
  interaction->ResetBit(kIAssumeFreeNucleon);

  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecNCEL::CalcChisq()
{
  double chisq = 0;
  chisq += TMath::Power(fMaTwkDial,   2.);
  chisq += TMath::Power(fEtaTwkDial,  2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNCEL::Init(void)
{
  AlgId id("genie::AhrensNCELPXSec","Default");

  AlgFactory * algf = AlgFactory::Instance();

  Algorithm * algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm * alg = algf->AdoptAlgorithm(id);
  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pDEBUG) << *fXSecModelConfig;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  this->SetMaPath ("Ma");
  this->SetEtaPath("Eta");

  fMaTwkDial   = 0.; 
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fEtaTwkDial  = 0.; 
  fEtaDef      = fXSecModelConfig->GetDouble(fEtaPath);
  fEtaCurr     = fEtaDef;

#ifdef _G_REWEIGHT_NCEL_DEBUG_
  fTestFile = new TFile("./ncel_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
