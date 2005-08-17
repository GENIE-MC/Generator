//____________________________________________________________________________
/*!

\class   genie::EventGenerator

\brief   Encapsulates a full ordered list of (is the aggregate of) concrete
         EventGeneratorI implementations that must act on the EventRecord
         to generate an event. Each of these implementations corresponds to
         a single processing step.

         Is a concrete implementation of the EventGeneratorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_H_
#define _EVENT_GENERATOR_H_

#include "EVGCore/EventGeneratorI.h"

namespace genie {

class EventGenerator : public EventGeneratorI {

public :

  EventGenerator();
  EventGenerator(const char * param_set);
  ~EventGenerator();

  //-- implement the original EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- implement the extensions to the EventRecordVisitorI interface

  const GVldContext &               ValidityContext  (void) const;
  const InteractionListGeneratorI * IntListGenerator (void) const;
  const XSecAlgorithmI *            CrossSectionAlg  (void) const;

  void  InstantiateValidityContext (void);

  //-- override part of the interface to allow the instantiation of its
  //   validity context (from an encoded XML config entry) during the
  //   configuration step

  virtual void Configure  (const Registry & config);
  virtual void Configure  (string param_set);

private:

  const EventRecordVisitorI * ProcessingStep(int istep) const;
};

}      // genie namespace

#endif // _EVENT_GENERATOR_H_
