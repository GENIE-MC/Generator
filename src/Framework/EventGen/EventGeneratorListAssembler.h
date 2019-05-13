//____________________________________________________________________________
/*!

\class   genie::EventGeneratorListAssembler

\brief   Assembles a list of all the EventGeneratorI subclasses that can be
         employed during a neutrino event generation job.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created January 25, 2004

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_LIST_ASSEMBLER_H_
#define _EVENT_GENERATOR_LIST_ASSEMBLER_H_

#include "Framework/Algorithm/Algorithm.h"

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
