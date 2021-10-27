//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cstring>

#include <RVersion.h>
#include <TClonesArray.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/HELepton/EventGen/HENuElGenerator.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif // __GENIE_PYTHIA6_ENABLED__

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::math;

//___________________________________________________________________________
HENuElGenerator::HENuElGenerator() :
EventRecordVisitorI("genie::HENuElGenerator")
{
  born = new Born();
}
//___________________________________________________________________________
HENuElGenerator::HENuElGenerator(string config) :
EventRecordVisitorI("genie::HENuElGenerator", config)
{

}
//___________________________________________________________________________
HENuElGenerator::~HENuElGenerator()
{

}
//___________________________________________________________________________
void HENuElGenerator::ProcessEventRecord(GHepRecord * event) const
{


  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo & proc_info   = interaction->ProcInfo();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu = event->Probe();
  GHepParticle * el = event->HitElectron();

  GHepParticle * target = event -> TargetNucleus();
  if(target) event->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  long double Ev = init_state.ProbeE(kRfLab); 
  long double mlout = interaction->FSPrimLepton()->Mass();
  long double mlout2  = mlout*mlout;
  
  long double s = 2 * kElectronMass * Ev + kElectronMass2;

  long double n1 = interaction->Kine().GetKV(kKVn1);
  long double n2 = interaction->Kine().GetKV(kKVn2);

  long double costh = n1;
  long double sinth = sqrtl(1-costh*costh);
  
  long double t = born->GetT(0.,kElectronMass,mlout,0.,s,n1);
  long double zeta  = born->GetReAlpha()/kPi*(2.0*logl(sqrtl(-t)/kElectronMass)-1.0);
  long double omx   = powl(n2, 1.0/zeta );
  long double s_r = s*( 1.-omx );

  // Boost velocity CM -> LAB
  long double EvCM = (s_r-kElectronMass2)/sqrtl(s_r)/2.;
  long double beta = (powl(Ev,2)-powl(EvCM,2))/(powl(Ev,2)+powl(EvCM,2));

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  long double Elpout = (s_r+mlout2)/sqrtl(s_r)/2.;
  long double Enuout = (s_r-mlout2)/sqrtl(s_r)/2.;
  LongLorentzVector p4_lpout( 0.,  Enuout*sinth,  Enuout*costh, Elpout );
  LongLorentzVector p4_nuout( 0., -Enuout*sinth, -Enuout*costh, Enuout );

  p4_lpout.BoostZ(beta);
  p4_nuout.BoostZ(beta);

  TLorentzVector p4lp_o( (double)p4_lpout.Px(), (double)p4_lpout.Py(), (double)p4_lpout.Pz(), (double)p4_lpout.E() );
  TLorentzVector p4nu_o( (double)p4_nuout.Px(), (double)p4_nuout.Py(), (double)p4_nuout.Pz(), (double)p4_nuout.E() );

  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2* kPi * rnd->RndLep().Rndm();
  p4lp_o.RotateZ(phi);
  p4nu_o.RotateZ(phi);

  //rotate from LAB=[0,0,Ev,Ev]->[px,py,pz,E]
  p4lp_o.RotateUz(unit_nu);
  p4nu_o.RotateUz(unit_nu);

  int pdgvout = 0;
  if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
  else if ( pdg::IsPositron(pdgl) ) pdgvout = kPdgNuE;
  else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
  else if ( pdg::IsAntiMuon(pdgl) ) pdgvout = kPdgNuMu;
  else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;
  else if ( pdg::IsAntiTau(pdgl)  ) pdgvout = kPdgNuTau;

  // Create a GHepParticle and add it to the event record
  event->AddParticle( pdgl,     kIStStableFinalState, 4, -1, -1, -1, p4lp_o,              *(nu->X4()) );
  event->AddParticle( pdgvout,  kIStStableFinalState, 4, -1, -1, -1, p4nu_o,              *(nu->X4()) );
  event->Summary()->KinePtr()->SetFSLeptonP4(p4lp_o);

}
//___________________________________________________________________________
void HENuElGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HENuElGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HENuElGenerator::LoadConfig(void)
{


}
//____________________________________________________________________________
