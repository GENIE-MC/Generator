//____________________________________________________________________________
/*!

\class   genie::EventRecord

\brief   Generated Event Record. It is a GHepRecord object that can accept /
         be visited by EventRecordVisitorI objects (event generation modules).
         All the other important container manipulation methods are defined
         at the base GHepRecord record.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _EVENT_RECORD_H_
#define _EVENT_RECORD_H_

#include <ostream>

#include "GHEP/GHepRecord.h"

using std::ostream;

namespace genie {

class EventRecordVisitorI;

class EventRecord : public GHepRecord {

public :

  EventRecord();
  EventRecord(int size);
  EventRecord(const EventRecord & record);
  ~EventRecord();

  void AcceptVisitor (EventRecordVisitorI * visitor);
  void Copy          (const EventRecord & record);
  void Print         (ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const EventRecord & event);

private:

ClassDef(EventRecord, 1)

};

}      // genie namespace

#endif // _EVENT_RECORD_H_
