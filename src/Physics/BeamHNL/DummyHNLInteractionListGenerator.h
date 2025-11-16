//____________________________________________________________________________
/*!

\class    genie::DummyHNLInteractionListGenerator

\brief    

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

	  Costas Andreopoulos <c.andreopoulos \at cern.ch>
	  University of Liverpool

\created  February 9th, 2022

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DUMMY_HNL_INTERACTION_LIST_GENERATOR_H_
#define _DUMMY_HNL_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class DummyHNLInteractionListGenerator : public InteractionListGeneratorI {

public :
  DummyHNLInteractionListGenerator();
  DummyHNLInteractionListGenerator(string config);
 ~DummyHNLInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _DUMMY_HNL_INTERACTION_LIST_GENERATOR_H_
