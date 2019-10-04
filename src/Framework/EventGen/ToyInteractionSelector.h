//____________________________________________________________________________
/*!

\class   genie::ToyInteractionSelector

\brief   Generates random interactions.

         This is a 'toy' InteractionSelectorI to be used in event generation
         testing / debugging. Not to be used in event generation for physics
         purposes.

         Is a concrete implementation of the InteractionSelectorI interface.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created December 05, 2004

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _TOY_INTERACTION_SELECTOR_H_
#define _TOY_INTERACTION_SELECTOR_H_

#include "Framework/EventGen/InteractionSelectorI.h"

namespace genie {

class ToyInteractionSelector : public InteractionSelectorI {

public :
  ToyInteractionSelector();
  ToyInteractionSelector(string config);
  ~ToyInteractionSelector();

  //! implement the InteractionSelectorI interface
  EventRecord * SelectInteraction
    (const InteractionGeneratorMap * igmp, const TLorentzVector & p4) const;
};

}      // genie namespace

#endif // _TOY_INTERACTION_SELECTOR_H_
