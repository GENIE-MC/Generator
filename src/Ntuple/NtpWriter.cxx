//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A utility class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TFolder.h>

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCPlainRecord.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCJobConfig.h"
#include "Ntuple/NtpMCJobEnv.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
NtpWriter::NtpWriter(NtpMCFormat_t fmt, Long_t runnu) :
fNtpFormat(fmt),
fRunNu(runnu),
fOutFile(0),
fOutTree(0),
fNtpMCPlainRecord(0),
fNtpMCEventRecord(0),
fNtpMCTreeHeader(0)
{
  LOG("NtpWriter", pNOTICE) << "Run number: " << runnu;
  LOG("NtpWriter", pNOTICE)
     << "Requested G/ROOT tree format: " << NtpMCFormat::AsString(fNtpFormat);
}
//____________________________________________________________________________
NtpWriter::~NtpWriter()
{
//  if(fNtpMCEvent) delete fNtpMCEvent;
//  if(fOutTree)    delete fOutTree;
//  if(fOutFile)    delete fOutFile;
//  if(fNtpMCTreeHeader) delete fNtpMCTreeHeader;
}
//____________________________________________________________________________
void NtpWriter::AddEventRecord(int ievent, const EventRecord * ev_rec)
{
//  if(fNtpMCEvent) delete fNtpMCEvent;

  LOG("NtpWriter", pINFO) << "Adding event " << ievent << " to output tree";

  switch (fNtpFormat) {
     case kNFPlainRecord:
          fNtpMCPlainRecord = new NtpMCPlainRecord();
          fNtpMCPlainRecord->Fill(ievent, ev_rec);
          break;
     case kNFEventRecord:
          fNtpMCEventRecord = new NtpMCEventRecord();
          fNtpMCEventRecord->Fill(ievent, ev_rec);
          break;
     default:
        break;
  }

  if(fOutTree) fOutTree->Fill();
  else {
    LOG("Ntp", pERROR) << "*** No TTree to add the input EventRecord!";
  }
}
//____________________________________________________________________________
void NtpWriter::Initialize(string filename_prefix)
{
  LOG("NtpWriter",pINFO) << "Initializing GENIE output MC tree";

  this->OpenFile(filename_prefix); // open ROOT file
  this->CreateTree();              // create output ROOT file

  TBranch * branch = this->CreateTreeBranch(); // create tree branch
  assert(branch);

  branch->SetAutoDelete(kFALSE);

  //-- create the tree header
  this->CreateTreeHeader();
  fNtpMCTreeHeader->Write();

  //-- save GENIE configuration for this MC Job
  NtpMCJobConfig configuration;
  configuration.Load()->Write();

  //-- take a snapshot of the user's environment
  NtpMCJobEnv environment;
  environment.TakeSnapshot()->Write();
}
//____________________________________________________________________________
void NtpWriter::OpenFile(string filename_prefix)
{
  if(fOutFile) delete fOutFile;

  // modify the filename to add the run number & the ntuple format
  ostringstream filename;
  filename << "GNtp-" << NtpMCFormat::FilenameTag(fNtpFormat)
           << "-" << fRunNu << ".root";

  LOG("NtpWriter", pINFO) << "Opening the output ROOT file: " << filename;
  fOutFile = new TFile(filename.str().c_str(),"RECREATE");
}
//____________________________________________________________________________
void NtpWriter::CreateTree(void)
{
  if(fOutTree) delete fOutTree;

  LOG("NtpWriter", pINFO) << "Creating the output GENIE/ROOT tree";

  ostringstream title;
  title << "GENIE MC Truth TTree"
              << ", Format: " << NtpMCFormat::AsString(fNtpFormat);

  fOutTree = new TTree("gtree",title.str().c_str());
  fOutTree->SetAutoSave(200000000);  // autosave when 0.2 Gbyte written
}
//____________________________________________________________________________
TBranch * NtpWriter::CreateTreeBranch(void)
{
  TBranch * branch = 0;

  switch (fNtpFormat) {
     case kNFPlainRecord:
        branch = this->CreatePRTreeBranch();
        break;
     case kNFEventRecord:
        branch = this->CreateERTreeBranch();
        break;
     default:
        LOG("NtpWriter", pERROR)
                    << "Unknown TTree format. Can not create TBranches";
        break;
  }
  return branch;
}
//____________________________________________________________________________
TBranch * NtpWriter::CreatePRTreeBranch(void)
{
  LOG("NtpWriter", pINFO) << "Creating a NtpMCPlainRecord TBranch";

  fNtpMCPlainRecord = 0;
  TTree::SetBranchStyle(2);

  TBranch * branch = fOutTree->Branch("gmcrec",
                      "genie::NtpMCPlainRecord", &fNtpMCPlainRecord, 32000,1);
  return branch;
}
//____________________________________________________________________________
TBranch * NtpWriter::CreateERTreeBranch(void)
{
  LOG("NtpWriter", pINFO) << "Creating a NtpMCEventRecord TBranch";

  fNtpMCEventRecord = 0;
  TTree::SetBranchStyle(1);

  TBranch * branch = fOutTree->Branch("gmcrec",
                    "genie::NtpMCEventRecord", &fNtpMCEventRecord, 32000, 1);
  return branch;
}
//____________________________________________________________________________
void NtpWriter::CreateTreeHeader(void)
{
  LOG("NtpWriter", pINFO) << "Creating the NtpMCTreeHeader";

  if(fNtpMCTreeHeader) delete fNtpMCTreeHeader;

  fNtpMCTreeHeader = new NtpMCTreeHeader;

  fNtpMCTreeHeader->format = fNtpFormat;
  fNtpMCTreeHeader->runnu  = fRunNu;

  LOG("NtpWriter", pINFO) << *fNtpMCTreeHeader;
}
//____________________________________________________________________________
void NtpWriter::Save(void)
{
  LOG("NtpWriter", pINFO) << "Saving the output tree";

  if(fOutFile) {

    fOutFile->Write();
    fOutFile->Close();
    delete fOutFile;
    fOutFile = 0;

  } else {
     LOG("NtpWriter", pERROR) << "No open ROOT file was found";
  }
}
//____________________________________________________________________________

