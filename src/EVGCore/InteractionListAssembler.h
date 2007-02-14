//____________________________________________________________________________
/*!

\class   genie::InteractionListAssembler

\brief   Assembles a list of all interactions that can be generated during a
         neutrino event generation job by querying each EventGeneratorI
         subclass employed in that job for its interaction list.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 16, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INTERACTION_LIST_ASSEMBLER_H_
#define _INTERACTION_LIST_ASSEMBLER_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class InteractionList;
class EventGeneratorList;
class InitialState;

class InteractionListAssembler : public Algorithm {

public :

  InteractionListAssembler();
  InteractionListAssembler(string config);
  ~InteractionListAssembler();

  void              SetGeneratorList        (EventGeneratorList * evglist);
  InteractionList * AssembleInteractionList (const InitialState & init) const;

private:

  EventGeneratorList * fEventGeneratorList;
};

}      // genie namespace

#endif // _INTERACTION_LIST_ASSEMBLER_H_
