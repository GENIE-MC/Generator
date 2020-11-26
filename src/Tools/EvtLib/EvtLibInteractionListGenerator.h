#ifndef _EVTLIB_INTERACTION_GENERATOR_H_
#define _EVTLIB_INTERACTION_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"

namespace genie {
namespace evtlib {

class EvtLibInteractionListGenerator : public InteractionListGeneratorI {

public :
  EvtLibInteractionListGenerator();
  EvtLibInteractionListGenerator(string config);
 ~EvtLibInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const override ;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);
};

} // evtlib namespace
} // genie namespace

#endif // _EVTLIB_INTERACTION_GENERATOR_H_
