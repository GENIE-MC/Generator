//____________________________________________________________________________
/*!

\class   genie::EventRecordVisitorI

\brief   Defines the EventRecordVisitorI interface.
         Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

*/
//____________________________________________________________________________

#ifndef _EVENT_RECORD_VISITOR_I_H_
#define _EVENT_RECORD_VISITOR_I_H_

#include "Algorithm/Algorithm.h"

namespace genie {

class GHepRecord;

class EventRecordVisitorI : public Algorithm {

public :

  //-- define the EventRecordVisitorI interface

  virtual void ProcessEventRecord(GHepRecord * event_rec) const = 0;

protected :

  EventRecordVisitorI();
  EventRecordVisitorI(const char * param_set);
  ~EventRecordVisitorI();
};

}      // genie namespace

#endif // _EVENT_RECORD_VISITOR_I_H_
