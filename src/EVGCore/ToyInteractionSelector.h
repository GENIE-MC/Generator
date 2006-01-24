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

#ifndef _TOY_INTERACTION_SELECTOR_H_
#define _TOY_INTERACTION_SELECTOR_H_

#include "EVGCore/InteractionSelectorI.h"

namespace genie {

class ToyInteractionSelector : public InteractionSelectorI {

public :

  ToyInteractionSelector();
  ToyInteractionSelector(string config);
  ~ToyInteractionSelector();

  //! implement the InteractionSelectorI interface
  EventRecord * SelectInteraction
          (const XSecAlgorithmMap * xsmp, const TLorentzVector & p4) const;
};

}      // genie namespace

#endif // _TOY_INTERACTION_SELECTOR_H_
