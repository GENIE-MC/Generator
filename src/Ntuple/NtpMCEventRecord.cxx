//____________________________________________________________________________
/*!

\class   genie::NtpMCEventRecord

\brief   MINOS-style ntuple record. Each such ntuple record holds a generated
         EventRecord object. Ntuples of this type are intended for feeding
         GENIE events into other applications (for example the GEANT4 based
         MC generation framework of an experiment) if no direct interface
         exists.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"

using std::endl;

using namespace genie;

ClassImp(NtpMCEventRecord)

//____________________________________________________________________________
EventRecord * genie::NtpMCEventRecord::fgevent = 0;
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
  this->Clear();

  this->event->Copy(*ev_rec);
  this->hdr.ievent = ievent;
}
//____________________________________________________________________________
void NtpMCEventRecord::Copy(const NtpMCEventRecord & ntpmcrec)
{
  this->Clear();

  this->event->Copy(*ntpmcrec.event);
  this->hdr.ievent = ntpmcrec.hdr.ievent;
}
//____________________________________________________________________________
void NtpMCEventRecord::Init(void)
{
  if(!fgevent) fgevent = new EventRecord;
  event = fgevent;
  hdr.ievent = 0;
}
//____________________________________________________________________________
void NtpMCEventRecord::Clear(void)
{
  event->Clear("C");
  hdr.ievent = 0;
}
//____________________________________________________________________________
