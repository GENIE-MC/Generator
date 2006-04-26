//____________________________________________________________________________
/*!

\class   genie::EGResponsibilityChain

\brief   A chain of EventGenerators.

         It implements a 'Chain of Responsibility' design pattern. \n
         The appropriate EventGenerator is selected based on the compatibility
         between its ModelValidityContext and the input Interaction.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GENERATOR_CHAIN_H_
#define _GENERATOR_CHAIN_H_

namespace genie {

class EventGeneratorI;
class EventGeneratorList;
class Interaction;

class EGResponsibilityChain {

public :

  EGResponsibilityChain();
  ~EGResponsibilityChain();

  void                    SetGeneratorList (EventGeneratorList * evglist);  
  const EventGeneratorI * FindGenerator    (const Interaction * interaction) const;

private:

  EventGeneratorList * fEventGeneratorList;
};

}      // genie namespace

#endif // _GENERATOR_CHAIN_H_
