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

#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCPlainRecord.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
NtpWriter::NtpWriter(NtpMCFormat_t fmt) :
fNtpFormat(fmt)
{
  this->Init();
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

  LOG("NtpWriter",pINFO) << "Adding event " << ievent << " to output tree";

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
  else
  {
    LOG("Ntp", pERROR) << "*** No TTree to add the input EventRecord!";
  }
}
//____________________________________________________________________________
void NtpWriter::InitTree(string filename)
{
  LOG("NtpWriter",pINFO) << "Initializing GENIE output MC tree";

  if(fOutFile) delete fOutFile;
  if(fOutTree) delete fOutTree;

  LOG("NtpWriter",pINFO) << "Creating ROOT file: " << filename;
  fOutFile = new TFile(filename.c_str(),"RECREATE");

  LOG("NtpWriter",pINFO)
             << "Creating a GENIE ROOT tree with format: "
                                       << NtpMCFormat::AsString(fNtpFormat);
  ostringstream title;
  title << "GENIE MC Truth TTree"
                       << ", Format: " << NtpMCFormat::AsString(fNtpFormat);

  fOutTree = new TTree("gtree",title.str().c_str());
  fOutTree->SetAutoSave(200000000);  // autosave when 0.2 Gbyte written

  LOG("NtpWriter",pINFO) << "Creating & saving the NtpMCTreeHeader";

  if(fNtpMCTreeHeader) delete fNtpMCTreeHeader;
  fNtpMCTreeHeader = new NtpMCTreeHeader;
  fNtpMCTreeHeader->Fill(fNtpFormat);
  LOG("NtpWriter",pINFO) << *fNtpMCTreeHeader;
  fNtpMCTreeHeader->Write();

  //TTree::SetBranchStyle(2);

  int split;
  int bufsize = 16000;
  TBranch *branch = 0;
  switch (fNtpFormat) {
     case kNFPlainRecord:
        LOG("NtpWriter",pINFO) << "Creating a NtpMCPlainRecord TBranch";

        fNtpMCPlainRecord = 0;
        split = 1;
        branch = fOutTree->Branch("gmcrec",
            "genie::NtpMCPlainRecord", &fNtpMCPlainRecord, bufsize, split);
        break;

     case kNFEventRecord:
        LOG("NtpWriter",pINFO) << "Creating a NtpMCEventRecord TBranch";
        fNtpMCEventRecord = 0;
        split = 99;
        branch = fOutTree->Branch("gmcrec",
            "genie::NtpMCEventRecord", &fNtpMCEventRecord, bufsize, split);
        break;

     default:
        LOG("NtpWriter", pERROR)
                          << "Unknown TTree format. Can not create TBranches";
        return;
        break;
  }
  branch->SetAutoDelete(kFALSE);
}
//____________________________________________________________________________
void NtpWriter::SaveTree(void)
{
  LOG("NtpWriter",pINFO) << "Saving the output tree";

  if(fOutFile) {

    fOutFile->Write();
    fOutFile->Close();
    delete fOutFile;
    fOutFile = 0;

  } else {
     LOG("NtpWriter",pERROR) << "No open ROOT file was found";
  }
}
//____________________________________________________________________________
void NtpWriter::Init(void)
{
  fOutFile          = 0;
  fOutTree          = 0;
  fNtpMCPlainRecord = 0;
  fNtpMCEventRecord = 0;
  fNtpMCTreeHeader  = 0;
}
//____________________________________________________________________________
