//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"

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

}
//___________________________________________________________________________
EventRecord::EventRecord(int size) :
GHepRecord(size)
{

}
//___________________________________________________________________________
EventRecord::EventRecord(const EventRecord & record) :
GHepRecord(record)
{

}
//___________________________________________________________________________
EventRecord::~EventRecord()
{

}
//___________________________________________________________________________
void EventRecord::AcceptVisitor(EventRecordVisitorI * visitor)
{
  visitor->ProcessEventRecord(this);
}
//___________________________________________________________________________
void EventRecord::Copy(const EventRecord & record)
{
  try {
    const GHepRecord & ghep = dynamic_cast<const GHepRecord &>(record);

    GHepRecord::Copy(ghep);

  } catch( std::bad_cast ) {
    LOG("EventRecord", pERROR) 
          << "Bad casting to 'const GHepRecord &'. Can not copy EventRecord";
  }
}
//___________________________________________________________________________
void EventRecord::Print(ostream & stream) const
{
  GHepRecord::Print(stream);
}
//___________________________________________________________________________
