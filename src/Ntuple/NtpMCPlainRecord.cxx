//____________________________________________________________________________
/*!

\class   genie::NtpMCPlainRecord

\brief   MINOS-style ntuple record. Each such ntuple record holds an MC event
         summary and a simplified version of the GHEP event record (as a
         TClones array of NtpMCGHepEntry objects). Ntuples of this type are
         intended for analysis in bare ROOT sessions.

\brief   MINOS-style Ntuple Class to a single generated MC event

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCPlainRecord.h"
#include "Ntuple/NtpMCGHepEntry.h"

using namespace genie;

ClassImp(NtpMCPlainRecord)

using std::endl;

//____________________________________________________________________________
TClonesArray * NtpMCPlainRecord::sghep = 0;
//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCPlainRecord & ntpp)
  {
     ntpp.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCPlainRecord::NtpMCPlainRecord():
NtpMCRecordI()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCPlainRecord::NtpMCPlainRecord(const NtpMCPlainRecord & ntpmcrec):
NtpMCRecordI()
{
  this->Copy(ntpmcrec);
}
//____________________________________________________________________________
NtpMCPlainRecord::~NtpMCPlainRecord()
{
  this->Clear();
}
//____________________________________________________________________________
void NtpMCPlainRecord::PrintToStream(ostream & stream) const
{
  NtpMCGHepEntry * ghep_entry = 0;
  TIter piter(ghep);
  while((ghep_entry = (NtpMCGHepEntry *) piter.Next())) stream << ghep_entry;
}
//____________________________________________________________________________
void NtpMCPlainRecord::Fill(unsigned int ievent, const EventRecord * ev_rec)
{
  TClonesArray & ghep_entries = *(this->ghep);

  int n = 0;
  GHepParticle * p = 0;
  TIter piter(ev_rec);

  while((p = (GHepParticle *) piter.Next())) {
    NtpMCGHepEntry * entry = new ( ghep_entries[n] ) NtpMCGHepEntry();
    entry->Copy(*p);
    entry->idx = n;
    n++;
  }

  // fill in the header
  this->hdr.ievent = ievent;
  this->nentries   = n;

  // fill in the summary
  this->mc.Copy(*ev_rec);
}
//____________________________________________________________________________
void NtpMCPlainRecord::Copy(const NtpMCPlainRecord & ntpmcrec)
{

}
//____________________________________________________________________________
void NtpMCPlainRecord::Init(void)
{
  if(!sghep) sghep = new TClonesArray("genie::NtpMCGHepEntry");
  sghep->SetOwner(true);

  this->ghep = sghep;
  this->mc.Init();
  this->nentries = 0;
  this->hdr.ievent = 0;
}
//____________________________________________________________________________
void NtpMCPlainRecord::Clear(void)
{
  sghep->Delete();

  this->ghep = 0;
  this->mc.Init();
  this->nentries = 0;
  this->hdr.ievent = 0;
}
//____________________________________________________________________________
