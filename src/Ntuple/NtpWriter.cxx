//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A simple class to facilitate creating the GENIE MC Ntuple from the
         output GENIE GHEP event records.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TTree.h>

#include "EventGeneration/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCRecord.h"

using namespace genie;

//____________________________________________________________________________
NtpWriter::NtpWriter()
{
  this->Init();
}
//____________________________________________________________________________
NtpWriter::~NtpWriter()
{
//  if(fNtpMCEvent) delete fNtpMCEvent;
//  if(fOutTree)    delete fOutTree;
//  if(fOutFile)    delete fOutFile;
}
//____________________________________________________________________________
void NtpWriter::AddEventRecord(int ievent, const EventRecord * ev_rec)
{
//  if(fNtpMCEvent) delete fNtpMCEvent;
  
  fNtpMCRecord = new NtpMCRecord();

  fNtpMCRecord->BuildRecord(ievent, ev_rec);

  if(fOutTree) fOutTree->Fill();
  else
  {
    LOG("Ntp", pERROR) << "*** No TTree to add the input EventRecord!";
  }  
}
//____________________________________________________________________________
void NtpWriter::InitTree(const char * filename)
{
  if(fOutFile) delete fOutFile;
  if(fOutTree) delete fOutTree;
  
  fOutFile = new TFile(filename,"RECREATE");  
  fOutTree = new TTree("genie","GENIE MC Truth TTree");

  fOutTree->SetAutoSave(200000000);  // autosave when 0.2 Gbyte written

  int split   = 1;
  int bufsize = 16000;
  
  fNtpMCRecord = 0;

  TTree::SetBranchStyle(2);

  TBranch *branch = fOutTree->Branch(
               "gmcrec", "genie::NtpMCRecord", &fNtpMCRecord, bufsize, split);

  branch->SetAutoDelete(kFALSE);
}  
//____________________________________________________________________________
void NtpWriter::SaveTree(void)
{
  if(fOutFile) {

    fOutFile->Write();
    fOutFile->Close();

    delete fOutFile;

    fOutFile = 0;    
  }
}
//____________________________________________________________________________
void NtpWriter::Init(void)
{
  fOutFile = 0;
  fOutTree = 0;
}
//____________________________________________________________________________
