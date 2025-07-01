//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "TBuffer.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

using std::endl;

using namespace genie;

ClassImp(NtpMCEventRecord)

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCEventRecord & ntpp)
  {
     ntpp.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCEventRecord::NtpMCEventRecord() :
NtpMCRecordI()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCEventRecord::NtpMCEventRecord(const NtpMCEventRecord & ntpmcrec) :
NtpMCRecordI()
{
  this->Copy(ntpmcrec);
}
//____________________________________________________________________________
NtpMCEventRecord::~NtpMCEventRecord()
{
  this->Clear();
}
//____________________________________________________________________________
void NtpMCEventRecord::PrintToStream(ostream & stream) const
{
  stream << this->hdr    << endl;
  stream << *this->event << endl;
}
//____________________________________________________________________________
void NtpMCEventRecord::Fill(unsigned int ievent, const EventRecord * ev_rec)
{
  this->event->Copy(*ev_rec);
  this->hdr.ievent = ievent;
}
//____________________________________________________________________________
void NtpMCEventRecord::Copy(const NtpMCEventRecord & ntpmcrec)
{
  this->event->Copy(*ntpmcrec.event);
  this->hdr.ievent = ntpmcrec.hdr.ievent;
}
//____________________________________________________________________________
void NtpMCEventRecord::Init(void)
{
  this->event      = new EventRecord;
  this->hdr.ievent = 0;
}
//____________________________________________________________________________
void NtpMCEventRecord::Clear(Option_t * /*opt*/)
{
  delete (this->event);
  this->event      = 0;
  this->hdr.ievent = 0;
}
//____________________________________________________________________________

// Extracted from the code ROOT generated
// added check-and-delete before streaming to
// avoid memory leak
void NtpMCEventRecord::Streamer(TBuffer &R__b) {
  // Stream an object of class genie::NtpMCEventRecord.

  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
    if (R__v) {
    }
    genie::NtpMCRecordI::Streamer(R__b);
    if (event) {
      delete event;
      event = nullptr;
    }
    R__b >> event;
    R__b.CheckByteCount(R__s, R__c, ::genie::NtpMCEventRecord::IsA());
  } else {
    R__c = R__b.WriteVersion(::genie::NtpMCEventRecord::IsA(), kTRUE);
    genie::NtpMCRecordI::Streamer(R__b);
    R__b << event;
    R__b.SetByteCount(R__c, kTRUE);
  }
}