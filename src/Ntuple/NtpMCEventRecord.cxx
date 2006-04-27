//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"

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
NtpMCEventRecord::NtpMCEventRecord()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCEventRecord::NtpMCEventRecord(const NtpMCEventRecord & ntpmcrec)
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
