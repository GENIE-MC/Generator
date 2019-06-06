//____________________________________________________________________________
/*!

\class    genie::SKInteractionListGenerator

\brief    Concrete implementations of the InteractionListGeneratorI interface.
          Creates a list of all the interactions that can be generated 
          by the single-Kaon generator.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 20, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SK_INTERACTION_LIST_GENERATOR_H_
#define _SK_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class SKInteractionListGenerator : public InteractionListGeneratorI {

public :

  SKInteractionListGenerator();
  SKInteractionListGenerator(string config);
 ~SKInteractionListGenerator();

  // Implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  bool fIsCC;
  bool fIsNC;
};

}      // genie namespace

#endif // _SK_INTERACTION_LIST_GENERATOR_H_
