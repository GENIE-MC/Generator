//____________________________________________________________________________
/*!

\class   genie::EventRecord

\brief   Generated Event Record. It is a GHepRecord object that can accept /
         be visited by EventRecordVisitorI objects (event generation modules).
         All the other important container manipulation methods are defined
         at the base GHepRecord record.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventRecord.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"

using namespace genie;

ClassImp(EventRecord)

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const EventRecord & event_record)
 {
   event_record.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
EventRecord::EventRecord() :
GHepRecord()
{
  fInteraction = 0;

  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::EventRecord(int size) :
GHepRecord(size)
{
  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::EventRecord(const EventRecord & record) :
GHepRecord(record)
{
  this->SetOwner(true);
}
//___________________________________________________________________________
EventRecord::~EventRecord()
{
  this->Delete();
}
//___________________________________________________________________________
void EventRecord::AcceptVisitor(EventRecordVisitorI * visitor)
{
  visitor->ProcessEventRecord(this);
}
//___________________________________________________________________________
void EventRecord::Print(ostream & stream) const
{
  GHepRecord::Print(stream);
}
//___________________________________________________________________________
