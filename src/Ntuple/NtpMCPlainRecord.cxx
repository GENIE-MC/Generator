//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
void NtpMCPlainRecord::Clear(Option_t * /*opt*/)
{
  sghep->Delete();

  this->ghep = 0;
  this->mc.Init();
  this->nentries = 0;
  this->hdr.ievent = 0;
}
//____________________________________________________________________________
