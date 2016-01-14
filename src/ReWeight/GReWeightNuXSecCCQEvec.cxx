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
 @ May 24, 2010 - CA
   First included in v2.7.1.

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
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"

//#define _G_REWEIGHT_CCQE_VEC_DEBUG_

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEvec::GReWeightNuXSecCCQEvec() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEvec::~GReWeightNuXSecCCQEvec()
{
#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_   
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEvec::IsHandled(GSyst_t syst)
{
   switch(syst) {
    case ( kXSecTwkDial_VecFFCCQEshape ) : 
       return true;
       break;
    default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_VecFFCCQEshape ) : 
       fFFTwkDial = twk_dial;
       break;
    default:
       return;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Reset(void)
{
  fFFTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEvec::CalcWeight(const genie::EventRecord & event) 
{
  bool tweaked = (TMath::Abs(fFFTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  Interaction * interaction = event.Summary();

  bool is_qe = interaction->ProcInfo().IsQuasiElastic();
  bool is_cc = interaction->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  bool charm = interaction->ExclTag().IsCharmEvent(); // skip CCQE charm
  if(charm) return 1.;

  int nupdg = event.Probe()->Pdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;
 
  //
  // Calculate weight
  // Input tweaking dial changes elastic nucleon form factors
  // (twk dial: 0 -> default/BBA, twk dial: 1 -> dipole).
  // Calculated weight includes `shape only effect in dsigma/dQ2
  // (normalized to constant integrated cross section)
  //

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double dial                = fFFTwkDial;
  double old_weight          = event.Weight();
  double def_xsec            = event.DiffXSec();
  double dpl_xsec            = fXSecModel_dpl->XSec(interaction, kPSQ2fE);
  double def_integrated_xsec = fXSecModel_bba->Integral(interaction);
  double dpl_integrated_xsec = fXSecModel_dpl->Integral(interaction);

  assert(def_integrated_xsec > 0.);
  assert(dpl_integrated_xsec > 0.);
//  if(def_integrated_xsec <= 0 || dpl_integrated_xsec <= 0) return 1.;

  double def_ratio = def_xsec / def_integrated_xsec;
  double dpl_ratio = dpl_xsec / dpl_integrated_xsec;

  assert(def_ratio > 0.);
//  if(def_ratio <= 0) return 1.;

  double weight = old_weight * (dial * dpl_ratio + (1-dial)*def_ratio) / def_ratio;

#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(
    E,Q2,weight,def_integrated_xsec,dpl_integrated_xsec,def_xsec,dpl_xsec);
#endif


  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEvec::CalcChisq()
{
  double chisq = TMath::Power(fFFTwkDial, 2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEvec::Init(void)
{
  AlgFactory * algf = AlgFactory::Instance();

  AlgId id0("genie::LwlynSmithQELCCPXSec","Default");
  Algorithm * alg0 = algf->AdoptAlgorithm(id0);
  fXSecModel_bba = dynamic_cast<XSecAlgorithmI*>(alg0);
  fXSecModel_bba->AdoptSubstructure();

  AlgId id1("genie::LwlynSmithQELCCPXSec","DipoleELFF");
  Algorithm * alg1 = algf->AdoptAlgorithm(id1);
  fXSecModel_dpl = dynamic_cast<XSecAlgorithmI*>(alg1);
  fXSecModel_dpl->AdoptSubstructure();

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fFFTwkDial = 0.;

#ifdef _G_REWEIGHT_CCQE_VEC_DEBUG_
  fTestFile = new TFile("./ccqevec_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght:sig0:sig:dsig0:dsig");
#endif

}
//_______________________________________________________________________________________

