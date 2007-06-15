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

#include <cassert>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TFolder.h>

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
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
  LOG("NtpWriter", pNOTICE) << "Adding event " << ievent << " to output tree";

  if(!ev_rec) {
    LOG("Ntp", pERROR) << "NULL input EventRecord!";
    return;
  }

  if(!fOutTree) {
    LOG("Ntp", pERROR) << "No open output TTree to add the input EventRecord!";
    return;
  }

  switch (fNtpFormat) {
     case kNFEventRecord:
          fNtpMCEventRecord = new NtpMCEventRecord();
          fNtpMCEventRecord->Fill(ievent, ev_rec);
          fOutTree->Fill();
          delete fNtpMCEventRecord;
          fNtpMCEventRecord = 0;
          break;
     default:
        break;
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
  filename << filename_prefix << "-" 
           << NtpMCFormat::FilenameTag(fNtpFormat)
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

