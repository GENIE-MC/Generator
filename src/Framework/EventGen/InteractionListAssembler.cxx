//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 16, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionListAssembler.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/EventGeneratorList.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

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
     LOG("IntLst", pERROR)
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

     LOG("IntLst", pINFO)
            << "\nQuerying EventGenerator: " << evgen->Id().Key()
                                            << " for its Interaction List";

     // ask the event generator to produce a list of all interaction it can
     // generate for the input initial state

     const InteractionListGeneratorI * intlistgen = evgen->IntListGenerator();

     InteractionList * intlist = intlistgen->CreateInteractionList(init_state); // should delete ?

     LOG("IntLst", pINFO) << "\nGot list:\n" << *intlist;

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
