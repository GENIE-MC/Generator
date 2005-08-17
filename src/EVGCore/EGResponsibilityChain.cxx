//____________________________________________________________________________
/*!

\class   genie::EGResponsibilityChain

\brief   A chain of EventGenerators.

         It implements a 'Chain of Responsibility' design pattern. \n
         The appropriate EventGenerator is selected based on the compatibility
         between its ModelValidityContext and the input Interaction.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include "EVGCore/EGResponsibilityChain.h"
#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
EGResponsibilityChain::EGResponsibilityChain()
{
  fEventGeneratorList = 0;
}
//___________________________________________________________________________
EGResponsibilityChain::~EGResponsibilityChain()
{

}
//___________________________________________________________________________
void EGResponsibilityChain::SetGeneratorList(EventGeneratorList * evglist)
{
  fEventGeneratorList = evglist;
}
//___________________________________________________________________________
const EventGeneratorI * EGResponsibilityChain::FindGenerator(
                                       const Interaction * interaction) const
{
  if(!fEventGeneratorList) {
     LOG("EVGChain", pERROR)
              << "\n*** NULL Generator List! "
                   << "Can not find an Event Generator for " << *interaction;
     return 0;
  }
  LOG("EVGChain", pINFO)
              << "Looking for an EventGenerator to generate the input event";

  // loop over all available event generators in this MC run
  EventGeneratorList::const_iterator evgliter;
  for(evgliter = fEventGeneratorList->begin();
                        evgliter != fEventGeneratorList->end(); ++evgliter) {

     const EventGeneratorI * evgen = *evgliter;

     LOG("EVGChain", pINFO)
           << "Trying........." << evgen->Name() << "/" << evgen->ParamSet();

     //-- Ask the validity context of this event generator & check whether
     //   it is the appropriate one to handle this event

     const GVldContext & vldc = evgen->ValidityContext();
     bool is_valid = vldc.IsValid(interaction);

     //-- Return the right event generator so that it can be asked to start
     //   generate the selected event
     if(is_valid) {
        LOG("EVGChain", pINFO)
           << " *** Selected: " << evgen->Name() << "/" << evgen->ParamSet();
        return evgen;
     }
  }
  LOG("EVGChain", pFATAL)
                 << "\n The requested event can not be generated - Aborting";
  LOG("EVGChain", pFATAL) << *interaction;

  assert(false);

  return 0;
}
//___________________________________________________________________________
