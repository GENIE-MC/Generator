//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
//#include "Framework/Numerical/RandomGen.h"
//#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
//#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Interaction/Interaction.h"
#include "Tools/VMC/EventLibraryInterface.h"

using namespace genie;
using namespace genie::vmc;

//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface() :
EventRecordVisitorI("genie::EventLibraryInterface")
{

}
//____________________________________________________________________________
EventLibraryInterface::EventLibraryInterface(string config) :
EventRecordVisitorI("genie::EventLibraryInterface",config)
{

}
//____________________________________________________________________________
EventLibraryInterface::~EventLibraryInterface()
{

}
//____________________________________________________________________________
void EventLibraryInterface::ProcessEventRecord(GHepRecord * event) const
{
// Get event summary constructed by GENIE
//
  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  TLorentzVector * probe_p4 = init_state.GetProbeP4(kRfLab);
  const TLorentzVector probe_v4(0.,0.,0.,0.);

  int probe_pdgc = init_state.ProbePdg();



  LOG("ELI", pINFO) << "Adding neutrino [pdgc = " << probe_pdgc << "]";

  event->AddParticle(probe_pdgc, kIStInitialState, -1,-1,-1,-1, *probe_p4, probe_v4);

  delete probe_p4;

  bool is_nucleus = init_state.Tgt().IsNucleus();
  if(is_nucleus) {
    int    tgt_A    = init_state.Tgt().A();
    int    tgt_Z    = init_state.Tgt().Z();
    int    tgt_pdgc = pdg::IonPdgCode(tgt_A, tgt_Z);
    double tgt_M    = PDGLibrary::Instance()->Find(tgt_pdgc)->Mass();

    LOG("ELI", pINFO)
         << "Adding nucleus [A = " << tgt_A << ", Z = " << tgt_Z
         << ", pdg = " << tgt_pdgc << "]";

    event->AddParticle(tgt_pdgc,kIStInitialState,-1,-1,-1,-1, 0,0,0,tgt_M, 0,0,0,0);
  }

  LOG("ELI", pNOTICE)
    << "Simulating NHL decay ";

//
// do your thing here...
//

}
//____________________________________________________________________________
void EventLibraryInterface::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void EventLibraryInterface::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void EventLibraryInterface::LoadConfig(void)
{
// here can read configuration from corresponding XML file, as usual

//
// ...
//

}
//___________________________________________________________________________
