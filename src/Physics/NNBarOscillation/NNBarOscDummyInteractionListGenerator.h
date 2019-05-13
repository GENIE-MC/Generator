//____________________________________________________________________________
/*!

\class    genie::NOscDummyInteractionListGenerator

\brief    

\author   Jeremy Hewes, Georgia Karagiorgi
          University of Manchester

\created  November, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NNBAR_OSC_DUMMY_INTERACTION_GENERATOR_H_
#define _NNBAR_OSC_DUMMY_INTERACTION_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {

class NNBarOscDummyInteractionListGenerator : public InteractionListGeneratorI {

public :
  NNBarOscDummyInteractionListGenerator();
  NNBarOscDummyInteractionListGenerator(string config);
 ~NNBarOscDummyInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _NNBAR_OSC_DUMMY_INTERACTION_GENERATOR_H_
