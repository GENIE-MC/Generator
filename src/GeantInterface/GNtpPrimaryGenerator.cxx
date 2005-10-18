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

#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <G4Event.hh>
#include <G4Types.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

#include "EVGCore/EventRecord.h"
#include "GeantInterface/GNtpPrimaryGenerator.h"
#include "Messenger/Messenger.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepOrder.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCFormat.h"

using std::vector;
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
void GNtpPrimaryGenerator::GeneratePrimaryVertex(G4Event* g4event)
{
// Get next GENIE's GHEP event record and convert it to GEANT's G4Event

  //-- read GENIE event from the input ROOT file

  EventRecord & gevt = *(this->ReadNextEvent());

  int np = gevt.GetEntries();
  if (np == 0) {
     LOG("GEANTi", pWARN)
                 << "Got empty GENIE event! Nothing to convert to G4Event";
     return;
  }

  //-- convert primary vertex position & time from GENIE/ROOT TLorentzVector
  //   to GEANT/CLHEP G4ThreeVector + G4double and create G4PrimaryVertex

  TLorentzVector * vtx = gevt.GetParticle(GHepOrder::ProbePosition())->V4();

  G4ThreeVector gvtx_pos  (vtx->X(), vtx->Y(), vtx->Z());
  G4double      gvtx_time (vtx->T());

  G4PrimaryVertex * g4vertex = new G4PrimaryVertex(gvtx_pos, gvtx_time);

  //-- create & store GEANT's primary particles

  vector<G4PrimaryParticle*> g4primaries(np);

  for(int ip = 0; ip < np; ip++) {

    // GEANT4's & GENIE's particles
    G4PrimaryParticle * g4particle = 0;
    GHepParticle *      particle   = gevt.GetParticle(ip);

    // convert GENIE's GHepParticle to GEANT's G4PrimaryParticle
    G4int    pdgc = particle->PdgCode();
    G4double px   = particle->Px();
    G4double py   = particle->Py();
    G4double pz   = particle->Pz();
    G4double m    = particle->P4()->M();

    g4particle = new G4PrimaryParticle(pdgc, px*GeV, py*GeV, pz*GeV);
    g4particle->SetMass(m*GeV);

    g4primaries[ip] = g4particle;
  }

  //-- create daughter-lists
  for(int ip = 0; ip < np; ip++) {
    GHepParticle * particle = gevt.GetParticle(ip);

    if(particle->HasDaughters()) {
       int child1 = particle->FirstDaughter();
       int child2 = particle->LastDaughter();
       for(int ic = child1; ic <= child2; ic++) {
         g4primaries[ip]->SetDaughter(g4primaries[ic]);
       }
    }
  }

  //-- add G4PrimaryParticles to G4PrimaryVertex
  //   (adding initial state particles only, as all daughter particles have
  //    been set in the daughter lists)

  for(int ip = 0; ip < np; ip++) {
    GHepParticle * particle = gevt.GetParticle(ip);
    GHepStatus_t ist = particle->Status();
    if(ist == kIStInitialState || ist == kIstNucleonTarget) {
         g4vertex->SetPrimary(g4primaries[ip]);
    }
  }

  //-- add G4PrimaryVertex to G4Event
  g4event->AddPrimaryVertex(g4vertex);

  // GEANT4 must have taken ownership of all created G4PrimaryParticles (?)
  // Will not clean the g4primaries vector
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

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
    G4Exception("GNtpPrimaryGenerator:: Can not open file.");
    return;
  }
  fFile = new TFile(filename.c_str(),"READ");

  //-- get tree header

  LOG("GEANTi", pDEBUG) << "Getting tree header";
  fTreeHdr = dynamic_cast <NtpMCTreeHeader *> ( fFile->Get("header") );
  LOG("GEANTi", pDEBUG) << "Input tree header: " << *fTreeHdr;

  if(fTreeHdr->format != kNFEventRecord) {
    G4Exception("GNtpPrimaryGenerator:: Incompatible event tree format.");
    return;
  }

  //-- get event TTree

  LOG("GEANTi", pDEBUG) << "Getting event tree";
  fTree = dynamic_cast <TTree *> (fFile->Get("gtree"));

  fNEntries = fTree->GetEntries();
  LOG("GEANTi", pINFO) << "Number of events to read: " << fNEntries;

  if(fNEntries <=0) {
    G4Exception("GNtpPrimaryGenerator:: Empty event tree.");
    return;
  }

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
    G4Exception("GNtpPrimaryGenerator:: End of GENIE's event tree.");
    return 0;
  }
  return 0;
}
//____________________________________________________________________________
