//____________________________________________________________________________
/*!

\program testReadGNtuple

\brief   A simple test program to illustrate how to use the GENIE ROOT output
         TTrees of both formats.

         To run: testReadGNtuple -f filename
         where the filename points to a ROOT file containing a GENIE output
         TTree

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 02, 2005
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
#include "Ntuple/NtpMCPlainRecord.h"
#include "Ntuple/NtpMCSummary.h"
#include "Ntuple/NtpMCGHepEntry.h"
#include "Messenger/Messenger.h"

using std::string;
using namespace genie;

bool   testCommandLineArgs(int argc);
string getRootFilename    (int argc, char ** argv);
bool   checkRootFilename  (string filename);
void   readNtpEventRecord (TTree * tree); // reads GNtpER.root files
void   readNtpPlainRecord (TTree * tree); // reads GNtpPR.root files
//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments and get the ROOT filename
  string filename = getRootFilename(argc,argv);

  if ( !testCommandLineArgs(argc)   ) return 1;
  if ( !checkRootFilename(filename) ) return 2;

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(filename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("Main", pINFO) << "Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;

  //-- read & print the TTree entries [depends on the format]
  LOG("Main", pINFO) << "This ntuple's format is : "
                             << NtpMCFormat::SuggestedUsage(format);
  switch(format) {
     case kNFPlainRecord:
       // simple ntuple for analysis in bare-ROOT sessions
       readNtpPlainRecord(tree);
       break;
     case kNFEventRecord:
       // ntuple for passing data to other applications (eg the MC
       // simulation framework of an experiment)
       readNtpEventRecord(tree);
       break;
     default:
       LOG("Main", pERROR) << "Cannot read this GENIE TTree format";
       return 3;
       break;
  }

  LOG("Main", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
bool testCommandLineArgs(int argc)
{
  if(argc!=3) {
   LOG("Main", pERROR) << "Not enough command line arguments";
   LOG("Main", pINFO)  << "Syntax: testReadGNtuple -f root_filename";
   return false;
  }
  return true;
}
//___________________________________________________________________
string getRootFilename(int argc, char ** argv)
{
  //-- Scan for filename from the command line argument (following -f)
  string filename="";
  for(int iarg = 0; iarg < argc-1; iarg++) {
   string argument(argv[iarg]);
   if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }
  return filename;
}
//___________________________________________________________________
bool checkRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("Main", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//___________________________________________________________________
void readNtpEventRecord(TTree * tree)
{
// This ntuple contains a single TBranch with NtpMCEventRecord
// objects in its leaves
// See NtpMCEventRecord for a description of ntuple contents
// This ntuple format is better used for moving GENIE data between
// applications (eg GENIE -> GEANT)

  if(!tree) return;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over TTree NtpMC records, get & print the events

  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("Main", pINFO) << rec_header;
    LOG("Main", pINFO) << event;

    mcrec->Clear();
  }
}
//___________________________________________________________________
void readNtpPlainRecord(TTree * tree)
{
// This ntuple contains a single TBranch with NtpMCPlainRecord
// objects in its leaves
// See NtpMCPlainRecord for a description of ntuple contents
// This ntuple format allows someone to analyze GENIE ntuples in
// bare ROOT sessions without loading the GENIE libraries.

  if(!tree) return;

  NtpMCPlainRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over TTree NtpMC records, get & print the events

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;    // record header
    NtpMCSummary   summary    = mcrec->mc;     // MC event summary
    TClonesArray & ghep       = *(mcrec->ghep);// Ntuple event record

    LOG("Main", pINFO) << rec_header;
    LOG("Main", pINFO) << "MC event summary: " << summary;

    // this loops over a simplified (ntuple, struct-like) version of
    // the EventRecord object which only contains public member data
    LOG("Main", pINFO) << "Ntuple-version of the EventRecord object:";
    LOG("Main", pINFO)
              << "|Idx | Ist |         PDG | mom | daughter  |"
                         <<"      Px |      Py |      Pz |       E |";

    for(int ientry = 0; ientry < ghep.GetEntries(); ientry++) {
      NtpMCGHepEntry * ghep_entry =
                   dynamic_cast<NtpMCGHepEntry *> (ghep[ientry]);
      LOG("Main", pINFO) << *ghep_entry;
    }
  }
}
//___________________________________________________________________

