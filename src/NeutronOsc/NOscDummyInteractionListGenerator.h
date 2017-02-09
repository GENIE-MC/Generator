//____________________________________________________________________________
/*!

\class    genie::NOscDummyInteractionListGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 10, 2011

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NOSC_DUMMY_INTERACTION_GENERATOR_H_
#define _NOSC_DUMMY_INTERACTION_GENERATOR_H_

#include "EVGCore/InteractionListGeneratorI.h"

namespace genie {

class NOscDummyInteractionListGenerator : public InteractionListGeneratorI {

public :
  NOscDummyInteractionListGenerator();
  NOscDummyInteractionListGenerator(string config);
 ~NOscDummyInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;
};

}      // genie namespace
#endif // _NOSC_DUMMY_INTERACTION_GENERATOR_H_
