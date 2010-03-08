//____________________________________________________________________________
/*!

\program gtestRwNuXSecHelperP2P

\brief   A simple program to test the GReWeightNuXSecHelper class.
	 The program will re-weight events generated/stored at some point
	 during the generator evolution to account for model changes in a
         subsequent improved / or tweaked / or bug-fixed version.
	 So reweighting happens between two {model/configuration} choices.

\syntax  gtestRwNuXSecHelperP2P -f filename
         where the filename points to a ROOT file with a GENIE event tree

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 10, 2007

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>

#include "Algorithm/AlgId.h"
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeightNuXSecHelper.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

string gOptInpFile; // input options (from command line arguments):

//___________________________________________________________________
int main(int argc, char ** argv)
{
  // Get the command line arguments
  GetCommandLineArgs(argc, argv);

  // Create a weight calculator
  rew::GReWeightNuXSecHelper wcalc;

  // Open the file and get the TTree & its header
  TFile file(gOptInpFile.c_str(),"READ");
  TTree * tree = dynamic_cast <TTree *> ( file.Get("gtree")  );
  NtpMCTreeHeader * thdr = 
       dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("test", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  assert(format == kNFGHEP); // only GHEP trees in this test

  // Set the branch address
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  // Loop over the event tree (GHEP records) and reweight events
  for(int i = 0; i< tree->GetEntries(); i++) 
  {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("test", pINFO) << rec_header;
    LOG("test", pINFO) << event;

    // reweight the event
    double wght = wcalc.NewWeight(event);

    LOG("test", pINFO)  
        << "Re-weighting: old wght. = " << event.Weight() 
        << ", new wght. = " << wght;

    mcrec->Clear();
  }

  file.Close();

  LOG("test", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // Get input ROOT file (containing a GENIE GHEP event tree)
  try {
    LOG("test", pINFO) << "Reading input filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("test", pFATAL)
             << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // Check input GENIE ROOT file
  bool ok = !(gSystem->AccessPathName(gOptInpFile.c_str()));
  if (!ok) {
    LOG("test", pFATAL)
      << "Input ROOT file [" << gOptInpFile << "] is not accessible";
    exit(2);
  }

  LOG("test", pNOTICE) 
      << "Input GENIE event file: " << gOptInpFile;
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("test", pNOTICE)
    << "\n\n" << "Syntax:" 
    << "\n gtestRwNuXSecHelperP2P -f ghep_file\n";
}
//___________________________________________________________________
