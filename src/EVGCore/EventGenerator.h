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

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_H_
#define _EVENT_GENERATOR_H_

#include <vector>

#include "EVGCore/EventGeneratorI.h"
#include "GHEP/GHepRecordHistory.h"

class TStopwatch;

using std::vector;

namespace genie {

class EventGenerator: public EventGeneratorI {

public :

  EventGenerator();
  EventGenerator(string config);
  ~EventGenerator();

  //-- implement the original EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- implement the extensions to the EventRecordVisitorI interface
  const GVldContext &               ValidityContext  (void) const;
  const InteractionListGeneratorI * IntListGenerator (void) const;
  const XSecAlgorithmI *            CrossSectionAlg  (void) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void Init                          (void);
  void LoadEVGModules                (void);
  void LoadInteractionListGenerator  (void);
  void LoadVldContext                (void);

  vector<const EventRecordVisitorI *> * fEVGModuleVec;
  vector<double> *                      fEVGTime;
  const XSecAlgorithmI *                fXSecModel;
  const InteractionListGeneratorI *     fIntListGen;
  GVldContext *                         fVldContext;
  TStopwatch *                          fWatch;
  mutable GHepRecordHistory             fRecHistory;

};

}      // genie namespace

#endif // _EVENT_GENERATOR_H_
