//____________________________________________________________________________
/*!

\class   genie::InteractionListAssembler

\brief   Assembles a list of all interactions that can be generated during a
         neutrino event generation job by querying each EventGeneratorI
         subclass employed in that job for its interaction list.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created May 16, 2005

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org     
*/
//____________________________________________________________________________

#ifndef _INTERACTION_LIST_ASSEMBLER_H_
#define _INTERACTION_LIST_ASSEMBLER_H_

#include "Framework/Algorithm/Algorithm.h"

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
