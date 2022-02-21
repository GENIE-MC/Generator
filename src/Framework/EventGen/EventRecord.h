//____________________________________________________________________________
/*!

\class   genie::EventRecord

\brief   Generated Event Record. It is a GHepRecord object that can accept /
         be visited by EventRecordVisitorI objects (event generation modules).
         All the other important container manipulation methods are defined
         at the base GHepRecord record.

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

\created October 1, 2004

\cpright Copyright (c) 2003-2022, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org         
*/
//____________________________________________________________________________

#ifndef _EVENT_RECORD_H_
#define _EVENT_RECORD_H_

#include <ostream>

#include "Framework/GHEP/GHepRecord.h"

using std::ostream;

namespace genie {

class EventRecord;
class EventRecordVisitorI;

ostream & operator<< (ostream& stream, const EventRecord & event);

class EventRecord : public GHepRecord {

public :
  using GHepRecord::Copy;  // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using GHepRecord::Print;

  EventRecord();
  EventRecord(int size);
  EventRecord(const EventRecord & record);
  ~EventRecord();

  void AcceptVisitor (EventRecordVisitorI * visitor);
  virtual void Copy          (const EventRecord & record);
  virtual void Print         (ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const EventRecord & event);

private:

ClassDef(EventRecord, 1)

};

}      // genie namespace

#endif // _EVENT_RECORD_H_
