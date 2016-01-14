//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2008 - CA
   Made event branch a priv data member and modified code accordingly.
   Added method to hand over the event tree so that additional user-defined
   branches (eg flux info) may be added.
 @ Mar 26, 2010 - CA
   Added CustomizeFilename() and CustomizeFilenamePrefix() to allow the use
   to customize either the entire output name or just the prefix before the
   run number.

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
fEventBranch(0),
fNtpMCEventRecord(0),
fNtpMCTreeHeader(0)
{
  LOG("Ntp", pNOTICE) << "Run number: " << runnu;
  LOG("Ntp", pNOTICE)
    << "Requested G/ROOT tree format: " << NtpMCFormat::AsString(fNtpFormat);

  this->SetDefaultFilename();
}
//____________________________________________________________________________
NtpWriter::~NtpWriter()
{

}
//____________________________________________________________________________
void NtpWriter::AddEventRecord(int ievent, const EventRecord * ev_rec)
{
  LOG("Ntp", pINFO) << "Adding event " << ievent << " to output tree";

  if(!ev_rec) {
    LOG("Ntp", pERROR) << "NULL input EventRecord!";
    return;
  }
  if(!fOutTree) {
    LOG("Ntp", pERROR) << "No open output TTree to add the input EventRecord!";
    return;
  }

  switch (fNtpFormat) {
     case kNFGHEP:
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
void NtpWriter::Initialize()
{
  LOG("Ntp",pINFO) << "Initializing GENIE output MC tree";

  this->OpenFile(fOutFilename); // open ROOT file
  this->CreateTree();           // create output tree

  //-- create the event branch
  this->CreateEventBranch(); 

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
void NtpWriter::CustomizeFilename(string filename)
{
 fOutFilename = filename;
}
//____________________________________________________________________________
void NtpWriter::CustomizeFilenamePrefix (string prefix)
{
  this->SetDefaultFilename(prefix);
}
//____________________________________________________________________________
void NtpWriter::SetDefaultFilename(string filename_prefix)
{
  ostringstream fnstr;
  fnstr << filename_prefix  << "." 
        << fRunNu << "."
        << NtpMCFormat::FilenameTag(fNtpFormat)
        << ".root";

  fOutFilename = fnstr.str();
}
//____________________________________________________________________________
void NtpWriter::OpenFile(string filename)
{
  if(fOutFile) delete fOutFile;

  LOG("Ntp", pINFO) 
      << "Opening the output ROOT file: " << filename;

  fOutFile = new TFile(filename.c_str(),"RECREATE");
}
//____________________________________________________________________________
void NtpWriter::CreateTree(void)
{
  if(fOutTree) delete fOutTree;

  LOG("Ntp", pINFO) << "Creating the output GENIE/ROOT tree";

  ostringstream title;
  title << "GENIE MC Truth TTree"
              << ", Format: " << NtpMCFormat::AsString(fNtpFormat);

  fOutTree = new TTree("gtree",title.str().c_str());
  fOutTree->SetAutoSave(200000000);  // autosave when 0.2 Gbyte written
}
//____________________________________________________________________________
void NtpWriter::CreateEventBranch(void)
{
  switch (fNtpFormat) {
     case kNFGHEP:
        this->CreateGHEPEventBranch();
        break;
     default:
        LOG("Ntp", pERROR)
           << "Unknown TTree format. Can not create TBranches";
        break;
  }
  assert(fEventBranch);
  fEventBranch->SetAutoDelete(kFALSE);
}
//____________________________________________________________________________
void NtpWriter::CreateGHEPEventBranch(void)
{
  LOG("Ntp", pINFO) << "Creating a NtpMCEventRecord TBranch";

  fNtpMCEventRecord = 0;
  TTree::SetBranchStyle(1);

  fEventBranch = fOutTree->Branch("gmcrec",
      "genie::NtpMCEventRecord", &fNtpMCEventRecord, 32000, 1);
}
//____________________________________________________________________________
void NtpWriter::CreateTreeHeader(void)
{
  LOG("Ntp", pINFO) << "Creating the NtpMCTreeHeader";

  if(fNtpMCTreeHeader) delete fNtpMCTreeHeader;

  fNtpMCTreeHeader = new NtpMCTreeHeader;

  fNtpMCTreeHeader->format = fNtpFormat;
  fNtpMCTreeHeader->runnu  = fRunNu;

  LOG("Ntp", pINFO) << *fNtpMCTreeHeader;
}
//____________________________________________________________________________
void NtpWriter::Save(void)
{
  LOG("Ntp", pINFO) << "Saving the output tree";

  if(fOutFile) {

    fOutFile->Write();
    fOutFile->Close();
    delete fOutFile;
    fOutFile = 0;

  } else {
     LOG("Ntp", pERROR) << "No open ROOT file was found";
  }
}
//____________________________________________________________________________

