//____________________________________________________________________________
/*!

\class   genie::EventGeneratorListAssembler

\brief   Assembles a list of all the EventGeneratorI subclasses that can be
         employed during a neutrino event generation job.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created January 25, 2004

*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_LIST_ASSEMBLER_H_
#define _EVENT_GENERATOR_LIST_ASSEMBLER_H_

#include "Algorithm/Algorithm.h"

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
