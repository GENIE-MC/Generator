//____________________________________________________________________________
/*!

\class   genie::InteractionListAssembler

\brief   Assembles a list of all interactions that can be generated during a
         neutrino event generation job by querying each EventGeneratorI
         subclass employed in that job for its interaction list.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 16, 2005

*/
//____________________________________________________________________________

#include "EVGCore/InteractionListAssembler.h"
#include "EVGCore/InteractionListGeneratorI.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
InteractionListAssembler::InteractionListAssembler() :
Algorithm("genie::InteractionListAssembler")
{
  fEventGeneratorList = 0;
}
//___________________________________________________________________________
InteractionListAssembler::InteractionListAssembler(string config) :
Algorithm("genie::InteractionListAssembler", config)
{
  fEventGeneratorList = 0;
}
//___________________________________________________________________________
InteractionListAssembler::~InteractionListAssembler()
{

}
//___________________________________________________________________________
void InteractionListAssembler::SetGeneratorList(EventGeneratorList * evglist)
{
  fEventGeneratorList = evglist;
}
//___________________________________________________________________________
InteractionList * InteractionListAssembler::AssembleInteractionList(
                                       const InitialState & init_state) const
{
  if(!fEventGeneratorList) {
     LOG("InteractionList", pERROR)
           << "\n*** NULL Generator List! "
             << "Can not assemble the Interaction List for \n" << init_state;
     return 0;
  }

  InteractionList * total_intlist = new InteractionList;

  EventGeneratorList::const_iterator evgliter; // event generator list iter
  InteractionList::iterator          intliter; // interaction list iter

  for(evgliter = fEventGeneratorList->begin();
                       evgliter != fEventGeneratorList->end(); ++evgliter) {

     const EventGeneratorI * evgen = *evgliter;

     LOG("InteractionList", pINFO)
            << "\nQuerying EventGenerator: " << evgen->Id().Key()
                                            << " for its Interaction List";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state

     const InteractionListGeneratorI * intlistgen = evgen->IntListGenerator();

     InteractionList * intlist = intlistgen->CreateInteractionList(init_state); // should delete ?

     LOG("InteractionList", pINFO) << "\nGot list:\n" << *intlist;

     // add them to the combined interaction list
     for(intliter = intlist->begin();
                                   intliter != intlist->end(); ++intliter) {
         Interaction * interaction = *intliter;
         total_intlist->push_back(interaction);

     } // loop over interaction that can be generated
  } // loop over event generators

  return total_intlist;
}
//___________________________________________________________________________
