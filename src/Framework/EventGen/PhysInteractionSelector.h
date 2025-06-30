//____________________________________________________________________________
/*!

\class   genie::PhysInteractionSelector

\brief   Selects interactions to be generated

         Is a concrete implementation of the InteractionSelectorI interface.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created January 25, 2005

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org 
*/
//____________________________________________________________________________

#ifndef _PHYS_INTERACTION_SELECTOR_H_
#define _PHYS_INTERACTION_SELECTOR_H_

#include "Framework/EventGen/InteractionSelectorI.h"

namespace genie {

class PhysInteractionSelector : public InteractionSelectorI {

public :
  PhysInteractionSelector();
  PhysInteractionSelector(string config);
  ~PhysInteractionSelector();

  //! implement the InteractionSelectorI interface
  EventRecord * SelectInteraction
     (const InteractionGeneratorMap * igmp, const TLorentzVector & p4) const;

  //! override the Algorithm::Configure methods to load configuration
  //! data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfigData (void);

  bool fUseSplines;
};

}      // genie namespace

#endif // _PHYS_INTERACTION_SELECTOR_H_
