//____________________________________________________________________________
/*!

\program gevdump

\brief   A GENIE utility printing-out GHEP event trees.

         Syntax:
           shell$ gevdump -f filename [-n nevents] [-o]

         [] denotes an optional argument

         -f Specifies a GENIE GHEP/ROOT event file.
         -n Specifies how many event to print (default: all)
         -o If set, instructs GENIE to print-out only the event whose number is
            specified via the `-n' option.

         Examples:

         1. Print out all events from /data/sample.ghep.root 
            shell$ gevdump -f /data/sample.ghep.root

         2. Print out the first 500 events from /data/sample.ghep.root 
            shell$ gevdump -f /data/sample.ghep.root -n 500

         3. Print out the event with id=178 from /data/sample.ghep.root 
            shell$ gevdump -f /data/sample.ghep.root -n 178 -o

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 02, 2005

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
bool CheckRootFilename  (string filename);
void Print1             (TTree * event_tree);
void PrintN             (TTree * event_tree);

int    gOptNEvt;
string gOptInpFilename;
bool   gOptPrint1;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  LOG("gevdump", pNOTICE) 
     << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  LOG("gevdump", pINFO) 
    << "This ntuple's format is : " << NtpMCFormat::AsString(format);

  if(gOptPrint1) Print1(tree);
  else           PrintN(tree);

  file.Close();

  LOG("gevdump", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void Print1(TTree * tree) 
{
  int nev = (int) tree->GetEntries();
  if(gOptNEvt < 0 || gOptNEvt >= nev) return;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  tree->GetEntry(gOptNEvt);

  NtpMCRecHeader rec_header = mcrec->hdr;
  EventRecord &  event      = *(mcrec->event);

  LOG("gevdump", pNOTICE) << rec_header;
  LOG("gevdump", pNOTICE) << event;

  mcrec->Clear();
}
//___________________________________________________________________
void PrintN(TTree * tree) 
{
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  //-- event loop
  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gevdump", pNOTICE) << rec_header;
    LOG("gevdump", pNOTICE) << event;

    mcrec->Clear();
  }
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevdump", pINFO) << "*** Parsing commad line arguments";

  // get GENIE event sample
  try {
    LOG("gevdump", pINFO) << "Reading event sample filename";
    gOptInpFilename = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevdump", pFATAL) 
        << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // number of events:
  try {    
    LOG("gevdump", pINFO) << "Reading number of events to analyze";
    gOptNEvt = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevdump", pINFO)
        << "Unspecified number of events to analyze - Use all";
      gOptNEvt = -1;
    }
  }

  gOptPrint1 = genie::utils::clap::CmdLineArgAsBool(argc,argv,'o');  
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevdump", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevdump -f sample.root [-n nev] [-o] \n";
}
//_________________________________________________________________________________
bool CheckRootFilename(string filename)
{
  if(filename.size() == 0) return false;
    
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gevdump", pERROR)  
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________

