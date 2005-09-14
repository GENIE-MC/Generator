//____________________________________________________________________________
/*!

\class   genie::geant::GNtpPrimaryGenerator

\brief   A GENIE/GEANT4 PrimaryGenerator. The NtpRdPrimaryGenerator feeds
         GENIE events into GEANT by simply reading them from GENIE's ER Tree.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created Sepember 12, 2005

*/
//____________________________________________________________________________

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include "EVGCore/EventRecord.h"
#include "GeantInterface/GNtpPrimaryGenerator.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"

using namespace genie;
using namespace genie::geant;

//____________________________________________________________________________
GNtpPrimaryGenerator::GNtpPrimaryGenerator()
{
  this->Initialize();
}
//____________________________________________________________________________
GNtpPrimaryGenerator::~GNtpPrimaryGenerator()
{
  this->CleanUp();
}
//____________________________________________________________________________
void GNtpPrimaryGenerator::GeneratePrimaryVertex(G4Event* geant_primary_event)
{
  //-- read GENIE event from the input ROOT file
  EventRecord * genie_event = this->ReadNextEvent();

  //-- convert GENIE's GHEP to GEANT's G4Event
}
//____________________________________________________________________________
void GNtpPrimaryGenerator::Initialize(void)
{
  fFile         = 0;
  fTree         = 0;
  fTreeHdr      = 0;
  fNtpRec       = 0;
  fNEntries     = 0;
  fCurrentEvent = 0;
}
//____________________________________________________________________________
void GNtpPrimaryGenerator::CleanUp(void)
{
  if(fFile) {
    fFile->Close();
    delete fFile;
    fFile = 0;
  }
  fTree         = 0;
  fTreeHdr      = 0;
  fNtpRec       = 0;
  fNEntries     = 0;
  fCurrentEvent = 0;
}
//____________________________________________________________________________
void GNtpPrimaryGenerator::ReadFromFile(string filename)
{
  //-- clean up
  this->CleanUp();

  //-- open ROOT file
  LOG("GEANTInterface", pINFO) << "Reading events from file: " << filename;

  fFile = new TFile(filename.c_str(),"READ");

  //-- get event TTree
  LOG("GEANTInterface", pDEBUG) << "Getting event tree";

  fTree = dynamic_cast <TTree *> (fFile->Get("gtree"));
  fNEntries = tree->GetEntries();
  LOG("GEANTInterface", pINFO) << "Number of events to read: " << fNEntries;

  //-- get tree header
  LOG("GEANTInterface", pDEBUG) << "Getting tree header";

  fTreeHdr = dynamic_cast <NtpMCTreeHeader *> ( fFile->Get("header") );
  LOG("GEANTInterface", pDEBUG) << "Input tree header: " << *fTreeHdr;

  //-- set the ntuple record
  LOG("GEANTInterface", pDEBUG) << "Setting NtpMCEventRecord branch addrress";

  tree->SetBranchAddress("gmcrec", &fNtpRec);
}
//____________________________________________________________________________
EventRecord * GNtpPrimaryGenerator::ReadNextEvent(void)
{
  if(fCurrentEvent < fNEntries) {

    LOG("GEANTInterface", pINFO) << "Getting tree entry: " << fCurrentEvent;
    fTree->GetEntry(fCurrentEvent);

    LOG("GEANTInterface", pINFO) << "Getting the stored GENIE event";
    NtpMCRecHeader rec_header = fNtpRec->hdr;
    EventRecord *  event      = fNtpRec->event;

    LOG("GEANTInterface", pINFO) << rec_header;
    LOG("GEANTInterface", pINFO) << *event;

    fCurrentEvent++;

    return event;

  } else {
    LOG("GEANTInterface", pWARN) << "No more events in the file!!!";
    return 0;
  }
  return 0;
}
//____________________________________________________________________________
