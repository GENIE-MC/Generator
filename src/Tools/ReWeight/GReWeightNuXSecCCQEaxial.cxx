//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Aaron Meyer <asmeyer \at uchicago.edu>
          University of Chicago/Fermilab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2015 - AM
   First version, based on GReWeightNuXSecCCQEvec.cxx

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>
#include <TFile.h>
#include <TNtupleD.h>

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
#include "Tools/ReWeight/GReWeightNuXSecCCQEaxial.h"
#include "Tools/ReWeight/GSystSet.h"
#include "Tools/ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecCCQEaxial::GReWeightNuXSecCCQEaxial() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecCCQEaxial::~GReWeightNuXSecCCQEaxial()
{
#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_   
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightNuXSecCCQEaxial::IsHandled(GSyst_t syst)
{
   switch(syst) {
    case ( kXSecTwkDial_AxFFCCQEshape ) : 
       return true;
       break;
    default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

  switch(syst) {
    case ( kXSecTwkDial_AxFFCCQEshape ) : 
       fFFTwkDial = twk_dial;
       break;
    default:
       return;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Reset(void)
{
  fFFTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEaxial::CalcWeight(const genie::EventRecord & event) 
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
  // (twk dial: 0 -> default/dipole, twk dial: 1 -> zexp).
  // Calculated weight includes `shape only effect in dsigma/dQ2
  // (normalized to constant integrated cross section)
  //

  interaction->KinePtr()->UseSelectedKinematics();
  interaction->SetBit(kIAssumeFreeNucleon);

  double dial                 = fFFTwkDial;
  double old_weight           = event.Weight();
  double def_xsec             = event.DiffXSec();
  double zexp_xsec            = fXSecModel_zexp->XSec(interaction, kPSQ2fE);
  //double def_integrated_xsec  = fXSecModel_dpl->Integral(interaction);
  //double zexp_integrated_xsec = fXSecModel_zexp->Integral(interaction);

  //assert(def_integrated_xsec > 0.);
  //assert(zexp_integrated_xsec > 0.);
//  if(def_integrated_xsec <= 0 || zexp_integrated_xsec <= 0) return 1.;

  //double def_ratio  = def_xsec  / def_integrated_xsec;
  //double zexp_ratio = zexp_xsec / zexp_integrated_xsec;
  double def_ratio  = def_xsec ;
  double zexp_ratio = zexp_xsec;

  assert(def_ratio > 0.);
//  if(def_ratio <= 0) return 1.;

  double weight = old_weight * (dial * zexp_ratio + (1-dial)*def_ratio) / def_ratio;

#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_
  double E  = interaction->InitState().ProbeE(kRfHitNucRest);
  double Q2 = interaction->Kine().Q2(true);
  fTestNtp->Fill(
    E,Q2,weight,def_integrated_xsec,zexp_integrated_xsec,def_xsec,zexp_xsec);
#endif


  return weight;
}
//_______________________________________________________________________________________
double GReWeightNuXSecCCQEaxial::CalcChisq()
{
  double chisq = TMath::Power(fFFTwkDial, 2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecCCQEaxial::Init(void)
{
  AlgFactory * algf = AlgFactory::Instance();

  AlgId id0("genie::LwlynSmithQELCCPXSec","Default");
  Algorithm * alg0 = algf->AdoptAlgorithm(id0);
  fXSecModel_dpl = dynamic_cast<XSecAlgorithmI*>(alg0);
  fXSecModel_dpl->AdoptSubstructure();

  AlgId id1("genie::LwlynSmithQELCCPXSec","ZExp");
  Algorithm * alg1 = algf->AdoptAlgorithm(id1);
  fXSecModel_zexp = dynamic_cast<XSecAlgorithmI*>(alg1);
  fXSecModel_zexp->AdoptSubstructure();

  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);

  fFFTwkDial = 0.;

#ifdef _G_REWEIGHT_CCQE_AXFF_DEBUG_
  fTestFile = new TFile("./ccqeaxil_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","E:Q2:wght:sig0:sig:dsig0:dsig");
#endif

}
//_______________________________________________________________________________________

