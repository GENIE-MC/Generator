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

#include <TLorentzVector.h>

#include "EVGCore/ToyInteractionSelector.h"
#include "EVGCore/EventRecord.h"
#include "EVGCore/InteractionList.h"
#include "EVGCore/InteractionFilter.h"
#include "EVGCore/XSecAlgorithmMap.h"
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

}
//___________________________________________________________________________
ToyInteractionSelector::ToyInteractionSelector(string config) :
InteractionSelectorI("genie::ToyInteractionSelector", config)
{

}
//___________________________________________________________________________
ToyInteractionSelector::~ToyInteractionSelector()
{

}
//___________________________________________________________________________
EventRecord * ToyInteractionSelector::SelectInteraction
          (const XSecAlgorithmMap * xscmap, const TLorentzVector & p4) const
{
  if(!xscmap) {
     LOG("InteractionSelector", pERROR)
               << "\n*** NULL XSecAlgorithmMap! Can't select interaction";
     return 0;
  }
  if(xscmap->size() <= 0) {
     LOG("InteractionSelector", pERROR)
              << "\n*** Empty XSecAlgorithmMap! Can't select interaction";
     return 0;
  }

  // select a random event generator
  RandomGen * rnd = RandomGen::Instance();

  const InteractionList & ilst = xscmap->GetInteractionList();

  unsigned int nint = ilst.size();
  unsigned int iint = (unsigned int) rnd->Random1().Integer(nint);

  Interaction * interaction = ilst[iint];

  // clone interaction
  Interaction * selected_interaction = new Interaction( *interaction );
  selected_interaction->GetInitialStatePtr()->SetProbeP4(p4);
  LOG("InteractionSelector", pINFO)
                   << "Interaction to generate: \n" << *selected_interaction;

  // bootstrap the event record
  EventRecord * evrec = new EventRecord;
  evrec->AttachInteraction(selected_interaction);

  return evrec;
}
//___________________________________________________________________________
