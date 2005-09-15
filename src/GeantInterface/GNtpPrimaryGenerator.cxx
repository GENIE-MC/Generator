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

#include <TFile.h>
#include <TTree.h>

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
void GNtpPrimaryGenerator::GeneratePrimaryVertex(G4Event* )
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
  LOG("GEANTi", pINFO) << "Reading events from file: " << filename;

  fFile = new TFile(filename.c_str(),"READ");

  //-- get event TTree
  LOG("GEANTi", pDEBUG) << "Getting event tree";

  fTree = dynamic_cast <TTree *> (fFile->Get("gtree"));
  fNEntries = fTree->GetEntries();
  LOG("GEANTi", pINFO) << "Number of events to read: " << fNEntries;

  //-- get tree header
  LOG("GEANTi", pDEBUG) << "Getting tree header";

  fTreeHdr = dynamic_cast <NtpMCTreeHeader *> ( fFile->Get("header") );
  LOG("GEANTi", pDEBUG) << "Input tree header: " << *fTreeHdr;

  //-- set the ntuple record
  LOG("GEANTi", pDEBUG) << "Setting NtpMCEventRecord branch addrress";

  fTree->SetBranchAddress("gmcrec", &fNtpRec);
}
//____________________________________________________________________________
EventRecord * GNtpPrimaryGenerator::ReadNextEvent(void)
{
  if(fCurrentEvent < fNEntries) {

    LOG("GEANTi", pINFO) << "Getting tree entry: " << fCurrentEvent;
    fTree->GetEntry(fCurrentEvent);

    LOG("GEANTi", pINFO) << "Getting the stored GENIE event";
    NtpMCRecHeader rec_header = fNtpRec->hdr;
    EventRecord *  event      = fNtpRec->event;

    LOG("GEANTi", pINFO) << rec_header;
    LOG("GEANTi", pINFO) << *event;

    fCurrentEvent++;

    return event;

  } else {
    LOG("GEANTi", pWARN) << "No more events in the file!!!";
    return 0;
  }
  return 0;
}
//____________________________________________________________________________
