//____________________________________________________________________________
/*!

\class   genie::EventGeneratorListAssembler

\brief   Assembles a list of all the EventGeneratorI subclasses that can be
         employed during a neutrino event generation job.

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created January 25, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org         
*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_LIST_ASSEMBLER_H_
#define _EVENT_GENERATOR_LIST_ASSEMBLER_H_

#include "Framework/Algorithm/Algorithm.h"
#include <string>
using std::string;

namespace genie {

class EventGeneratorList;
class EventGeneratorI;

class EventGeneratorListAssembler : public Algorithm {

public :

  EventGeneratorListAssembler();
  EventGeneratorListAssembler(string config);
  ~EventGeneratorListAssembler();

  EventGeneratorList * AssembleGeneratorList();

private:

  const EventGeneratorI * LoadGenerator(int ip);
};

}      // genie namespace

#endif // _EVENT_GENERATOR_LIST_ASSEMBLER_H_
