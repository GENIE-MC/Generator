//____________________________________________________________________________
/*!

\class    genie::DummyInteractionListGenerator

\brief

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  November 10, 2011

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _DUMMY_INTERACTION_GENERATOR_H_
#define _DUMMY_INTERACTION_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class DummyInteractionListGenerator : public InteractionListGeneratorI {

public :
  DummyInteractionListGenerator();
  DummyInteractionListGenerator(string config);
 ~DummyInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _DUMMY_INTERACTION_GENERATOR_H_
