//____________________________________________________________________________
/*!

\class   genie::ToyInteractionSelector

\brief   Generates random interactions.

         This is a 'toy' InteractionSelectorI to be used in event generation
         testing / debugging. Not to be used in event generation for physics
         purposes.

         Is a concrete implementation of the InteractionSelectorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 05, 2004

*/
//____________________________________________________________________________

#include "EVGCore/ToyInteractionSelector.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionFilter.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
ToyInteractionSelector::ToyInteractionSelector() :
InteractionSelectorI("genie::ToyInteractionSelector")
{
  fEventGeneratorList = 0;
  fInteractionFilter  = 0;
}
//___________________________________________________________________________
ToyInteractionSelector::ToyInteractionSelector(string config) :
InteractionSelectorI("genie::ToyInteractionSelector", config)
{
  fEventGeneratorList = 0;
  fInteractionFilter  = 0;
}
//___________________________________________________________________________
ToyInteractionSelector::~ToyInteractionSelector()
{

}
//___________________________________________________________________________
void ToyInteractionSelector::SetGeneratorList(
                                          const EventGeneratorList * evglist)
{
  fEventGeneratorList = evglist;
}
//___________________________________________________________________________
void ToyInteractionSelector::SetInteractionFilter(
                                            const InteractionFilter * filter)
{
  fInteractionFilter = filter;
}
//___________________________________________________________________________
Interaction * ToyInteractionSelector::SelectInteraction(
                                       const InitialState & init_state) const
{
  if(!fEventGeneratorList) {
     LOG("InteractionSelector", pERROR)
               << "\n*** NULL Generator List! "
                         << "Can not select interaction for " << init_state;
     return 0;
  }
  if(fEventGeneratorList->size() <= 0) {
     LOG("InteractionSelector", pERROR)
               << "\n*** Empty Generator List! "
                         << "Can not select interaction for " << init_state;
     return 0;
  }

  // select a random event generator
  RandomGen * rnd = RandomGen::Instance();

  unsigned int ngen = fEventGeneratorList->size();
  int          igen = rnd->Random1().Integer(ngen);

  const EventGeneratorI * evgen = (*fEventGeneratorList)[igen];

  // ask the event generator to produce a list of all interaction it can
  // generate for the input initial state

  const InteractionListGeneratorI * intlistgen = evgen->IntListGenerator();

  InteractionList * intlist = intlistgen->CreateInteractionList(init_state);

  // select a random interaction from the interaction list

  unsigned int nint = intlist->size();
  int          iint = rnd->Random1().Integer(nint);

  Interaction * interaction = (*intlist)[iint];

  // clone, print and return

  Interaction * selected_interaction = new Interaction( *interaction );

  LOG("InteractionSelector", pINFO)
                   << "Interaction to generate: \n" << *selected_interaction;

  delete intlist;
  return selected_interaction;
}
//___________________________________________________________________________
