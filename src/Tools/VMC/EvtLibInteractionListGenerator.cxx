#include "Tools/VMC/EvtLibInteractionListGenerator.h"

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::vmc;

//___________________________________________________________________________
EvtLibInteractionListGenerator::EvtLibInteractionListGenerator() :
InteractionListGeneratorI("genie::vmc::EvtLibInteractionListGenerator")
{

}
//___________________________________________________________________________
EvtLibInteractionListGenerator::EvtLibInteractionListGenerator(string config):
InteractionListGeneratorI("genie::vmc::EvtLibInteractionListGenerator",  config)
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

  // REVIEW these are all these arbitrary choices, but seems to be required to
  // be able to succesfully pick a spline. Perhaps it would be better to add
  // new enum entries for "unspecified because it came from an external
  // library"?

  ProcessInfo proc_info_cc(kScQuasiElastic, kIntWeakCC);
  Interaction* interaction_cc = new Interaction(init_state, proc_info_cc);
  interaction_cc->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgNeutron);
  intlist->push_back(interaction_cc);

  ProcessInfo proc_info_nc(kScQuasiElastic, kIntWeakNC);
  Interaction* interaction_nc = new Interaction(init_state, proc_info_nc);
  interaction_nc->InitStatePtr()->TgtPtr()->SetHitNucPdg(kPdgProton);
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
