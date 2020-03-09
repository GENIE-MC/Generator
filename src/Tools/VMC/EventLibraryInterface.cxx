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
//#include "Framework/ParticleData/PDGUtils.h"
//#include "Framework/ParticleData/PDGLibrary.h"
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



  LOG("ELI", pINFO) << "Adding neutrino [pdgc = " << pdgc << "]";

  event->AddParticle(pdgc,kIStInitialState, -1,-1,-1,-1, *p4, v4);

  delete p4;

  bool is_nucleus = init_state.Tgt().IsNucleus();
  if(is_nucleus) {
    int    tgt_A    = init_state.Tgt().A();
    int    tgt_Z    = init_state.Tgt().Z();
    int    tgt_pdgc = pdg::IonPdgCode(A, Z);
    double tgt_M    = PDGLibrary::Instance()->Find(pdgc)->Mass();

    LOG("ELI", pINFO)
         << "Adding nucleus [A = " << A << ", Z = " << Z
         << ", pdg = " << pdgc << "]";

  event->AddParticle(pdgc,kIStInitialState,-1,-1,-1,-1, 0,0,0,M, 0,0,0,0);


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
