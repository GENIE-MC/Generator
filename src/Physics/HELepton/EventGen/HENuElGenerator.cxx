//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
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
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu = event->Probe();

  GHepParticle * target = event -> TargetNucleus();
  if(target) event->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  bool isCC = proc_info.IsWeakCC();

  long double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton
  long double mlin  = kElectronMass;                       //mass of incoming charged lepton

  long double Enuin = init_state.ProbeE(kRfLab);
  long double s = born->GetS(mlin,Enuin);

  long double n1 = interaction->Kine().GetKV(kKVn1);
  long double n2 = interaction->Kine().GetKV(kKVn2);

  long double costhCM = n1;
  long double sinthCM = sqrtl(1-costhCM*costhCM);
  
  long double t  = born->GetT( mlin, mlout, s, n1 );
  long double zeta  = born->GetReAlpha()/kPi*(2.0*logl(sqrtl(-t)/kElectronMass)-1.0);
  long double omx   = powl(n2, 1.0/zeta );
  long double s_r = s*( 1.-omx );

  // Boost velocity CM -> LAB
  long double EnuinCM = (s_r-mlin*mlin)/sqrtl(s_r)/2.;
  long double beta = (powl(Enuin,2)-powl(EnuinCM,2))/(powl(Enuin,2)+powl(EnuinCM,2));

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  long double ElpoutCM = (s_r+mlout*mlout)/sqrtl(s_r)/2.;
  long double EnuoutCM = (s_r-mlout*mlout)/sqrtl(s_r)/2.;
  LongLorentzVector p4_lpout( 0.,  EnuoutCM*sinthCM,  EnuoutCM*costhCM, ElpoutCM );
  LongLorentzVector p4_nuout( 0., -EnuoutCM*sinthCM, -EnuoutCM*costhCM, EnuoutCM );

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
  if (isCC) pdgvout = kPdgNuE;
  else      pdgvout = nu->Pdg();

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
