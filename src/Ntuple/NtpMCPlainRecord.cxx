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
TClonesArray * genie::NtpMCPlainRecord::fgGHep = 0;
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
  this->Clear();

  // copy all GHepParticles

  TClonesArray & ghep_entries = *ghep;

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
  hdr.ievent   = ievent;
  nentries = n;

  // fill in the summary
  mc.Copy(*ev_rec);
}
//____________________________________________________________________________
void NtpMCPlainRecord::Copy(const NtpMCPlainRecord & ntpmcrec)
{

}
//____________________________________________________________________________
void NtpMCPlainRecord::Init(void)
{
  if(!fgGHep) fgGHep = new TClonesArray("genie::NtpMCGHepEntry");
  ghep = fgGHep;
  hdr.ievent = 0;
}
//____________________________________________________________________________
void NtpMCPlainRecord::Clear(void)
{
  ghep->Clear("C");
  nentries = 0;
  hdr.ievent = 0;
}
//____________________________________________________________________________
