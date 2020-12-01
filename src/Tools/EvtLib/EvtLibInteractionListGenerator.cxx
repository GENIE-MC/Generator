#include "Tools/EvtLib/EvtLibInteractionListGenerator.h"

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::evtlib;

//___________________________________________________________________________
EvtLibInteractionListGenerator::EvtLibInteractionListGenerator() :
InteractionListGeneratorI("genie::evtlib::EvtLibInteractionListGenerator")
{

}
//___________________________________________________________________________
EvtLibInteractionListGenerator::EvtLibInteractionListGenerator(string config):
InteractionListGeneratorI("genie::evtlib::EvtLibInteractionListGenerator",  config)
{

}
//___________________________________________________________________________
EvtLibInteractionListGenerator::~EvtLibInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * EvtLibInteractionListGenerator::CreateInteractionList(
   const InitialState & init_state) const
{
  InteractionList * intlist = new InteractionList;

  // Struck nucleon is arbitrary but seems to be required

  ProcessInfo proc_info_cc(kScUnknown, kIntWeakCC);
  Interaction* interaction_cc = new Interaction(init_state, proc_info_cc);
  intlist->push_back(interaction_cc);

  ProcessInfo proc_info_nc(kScUnknown, kIntWeakNC);
  Interaction* interaction_nc = new Interaction(init_state, proc_info_nc);
  intlist->push_back(interaction_nc);

  return intlist;
}

//____________________________________________________________________________
void EvtLibInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
}
//____________________________________________________________________________
void EvtLibInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
}
