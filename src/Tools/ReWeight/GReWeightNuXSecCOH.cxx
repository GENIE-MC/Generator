//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
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

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Tools/ReWeight/GReWeightNuXSecCOH.h"
#include "Tools/ReWeight/GSystSet.h"
#include "Tools/ReWeight/GSystUncertainty.h"
#include "Framework/Registry/Registry.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCOH::GReWeightNuXSecCOH() :
GReWeightModel("CCCoh"),
fManualModelName(),
fManualModelType()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCOH::GReWeightNuXSecCOH(std::string model, std::string type) :
GReWeightModel("CCCoh"),
fManualModelName(model),
fManualModelType(type)
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCOH::~GReWeightNuXSecCOH()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecCOH::IsHandled(GSyst_t syst)
{
   bool handle;

   switch(syst) {
     case ( kXSecTwkDial_MaCOHpi ) : 
     case ( kXSecTwkDial_R0COHpi ) : 
       handle = true;  
       break;
     default:
       handle = false;
       break;
   }
   return handle;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::SetSystematic(GSyst_t syst, double twk_dial)
{
   switch(syst) {
     case ( kXSecTwkDial_MaCOHpi ) : 
       fMaTwkDial = twk_dial;
       break;
     case ( kXSecTwkDial_R0COHpi ) : 
       fR0TwkDial = twk_dial;
       break;
     default:
       break;
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Reset(void)
{
  fMaTwkDial   = 0.; 
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.; 
  fR0Curr      = fR0Def;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Reconfigure(void)
{
  GSystUncertainty * fracerr = GSystUncertainty::Instance();

  double fracerr_ma = fracerr->OneSigmaErr(kXSecTwkDial_MaCOHpi);
  double fracerr_r0 = fracerr->OneSigmaErr(kXSecTwkDial_R0COHpi);

  fMaCurr = fMaDef * (1. + fMaTwkDial * fracerr_ma);
  fR0Curr = fR0Def * (1. + fR0TwkDial * fracerr_r0);

  fMaCurr = TMath::Max(0., fMaCurr  );
  fR0Curr = TMath::Max(0., fR0Curr  );

  Registry & r = *fXSecModelConfig;

  r.Set(fMaPath, fMaCurr); 
  r.Set(fR0Path, fR0Curr); 

  fXSecModel->Configure(r);

//LOG("ReW", pDEBUG) << *fXSecModel;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCOH::CalcWeight(const genie::EventRecord & event) 
{
  Interaction * interaction = event.Summary();

  bool is_coh = interaction->ProcInfo().IsCoherent();
  if(!is_coh) return 1.;

  bool tweaked = 
      (TMath::Abs(fMaTwkDial) > controls::kASmallNum) ||
      (TMath::Abs(fR0TwkDial) > controls::kASmallNum);
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

  interaction->KinePtr()->UseSelectedKinematics();
  
  const KinePhaseSpace_t phase_space = kPSxyfE;

  double old_xsec   = event.DiffXSec();
  if (!fUseOldWeightFromFile || fNWeightChecksDone < fNWeightChecksToDo) {
    double calc_old_xsec = fXSecModelDef->XSec(interaction, phase_space);
    if (fNWeightChecksDone < fNWeightChecksToDo) {
      if (std::abs(calc_old_xsec - old_xsec)/old_xsec > controls::kASmallNum) {
        LOG("ReW",pWARN) << "Warning - default dxsec does not match dxsec saved in tree. Does the config match?";
      }
      fNWeightChecksDone++;
    }
    if(!fUseOldWeightFromFile) {
      old_xsec = calc_old_xsec;
    }
  }
  
  double old_weight = event.Weight();
  double new_xsec   = fXSecModel->XSec(interaction, phase_space);
  double new_weight = old_weight * (new_xsec/old_xsec);

  interaction->KinePtr()->ClearRunningValues();

  return new_weight;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCOH::Init(void)
{
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * gpl = conf_pool->GlobalParameterList();
  RgAlg xsec_alg = gpl->GetAlg("XSecModel@genie::EventGenerator/COH-CC");
  
  AlgId id(xsec_alg);
  
  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory * algf = AlgFactory::Instance();
  Algorithm * alg_def = algf->AdoptAlgorithm(id);
  Algorithm * alg_twk = algf->AdoptAlgorithm(twk_id);

  fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg_twk);
  fXSecModel->AdoptSubstructure();
  
  fXSecModelDef = dynamic_cast<XSecAlgorithmI*>(alg_def);
  fXSecModelDef->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
//LOG("ReW", pNOTICE) << *fXSecModelConfig;

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  this->SetMaPath("Ma");
  this->SetR0Path("Ro");

  fMaTwkDial   = 0.; 
  fMaDef       = fXSecModelConfig->GetDouble(fMaPath);
  fMaCurr      = fMaDef;
  fR0TwkDial   = 0.; 
  fR0Def       = fXSecModelConfig->GetDouble(fR0Path);
  fR0Curr      = fR0Def;
}
//_______________________________________________________________________________________


