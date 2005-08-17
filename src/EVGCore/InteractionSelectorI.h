//____________________________________________________________________________
/*!

\class   genie::InteractionSelectorI

\brief   Defines the InteractionSelectorI interface to be implemented by
         algorithms selecting interactions to be generated.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 05, 2004

*/
//____________________________________________________________________________

#ifndef _INTERACTION_SELECTOR_I_H_
#define _INTERACTION_SELECTOR_I_H_

#include "Algorithm/Algorithm.h"
#include "Interaction/InitialState.h"

namespace genie {

class Interaction;
class InteractionFilter;
class EventGeneratorList;

class InteractionSelectorI : public Algorithm {

public :

  virtual ~InteractionSelectorI();

  //!  Define the InteractionSelectorI interface

  virtual void          SetInteractionFilter (const InteractionFilter * filt) = 0;
  virtual void          SetGeneratorList     (const EventGeneratorList * egl) = 0;
  virtual Interaction * SelectInteraction    (const InitialState & ist) const = 0;

protected:

  InteractionSelectorI();
  InteractionSelectorI(const char * param_set);
};

}      // genie namespace

#endif // _INTERACTION_SELECTOR_I_H_
