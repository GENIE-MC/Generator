//____________________________________________________________________________
/*!

\class   genie::PhysInteractionSelector

\brief   Selects interactions to be generated

         Is a concrete implementation of the InteractionSelectorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created January 25, 2005

*/
//____________________________________________________________________________

#ifndef _PHYS_INTERACTION_SELECTOR_H_
#define _PHYS_INTERACTION_SELECTOR_H_

#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/InteractionSelectorI.h"

namespace genie {

class InteractionFilter;
class EventGeneratorList;

class PhysInteractionSelector : public InteractionSelectorI {

public :

  PhysInteractionSelector();
  PhysInteractionSelector(string config);
  ~PhysInteractionSelector();

  //! implement the InteractionSelectorI interface

  void          SetInteractionFilter (const InteractionFilter * filt);
  void          SetGeneratorList     (const EventGeneratorList * egl);
  Interaction * SelectInteraction    (const InitialState & init_state) const;

private:

  bool  UseStoredCrossSections (void) const;

  const InteractionFilter *  fInteractionFilter;
  const EventGeneratorList * fEventGeneratorList;
};

}      // genie namespace

#endif // _PHYS_INTERACTION_SELECTOR_H_
