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
#include "Framework/Messenger/Messenger.h"
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
