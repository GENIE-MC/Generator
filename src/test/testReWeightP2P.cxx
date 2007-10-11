//____________________________________________________________________________
/*!

\program gtestReWeightP2P

\brief   A simple test program to illustrate how to use the GENIE event
         reweighting package.
	 The program will re-weight events generated/stored at some point
	 during the generator evolution to account for model changes in a
         subsequent improved / or tweaked / or bug-fixed version.
	 So reweighting happens between two {model/configuration} choices.
         Only cross section model reweighting is included.

\syntax  gtestReWeightP2P -f filename
         where the filename points to a ROOT file with a GENIE event tree

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 10, 2007

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
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
#include "ReWeight/ReWeightCrossSection.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

string gOptInpFile; // input options (from command line arguments):

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  //-- create a weight calculator
  ReWeightCrossSection wcalc;

  //-- open the file and get the TTree & its header
  TFile file(gOptInpFile.c_str(),"READ");
  TTree * tree = dynamic_cast <TTree *> ( file.Get("gtree")  );
  NtpMCTreeHeader * thdr = 
       dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("test", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord); // only ER trees in this test

  //-- set the branch address
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over the event tree (GHEP records) and reweight events

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
  //get input ROOT file (containing a GENIE ER/GHEP ntuple)
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

  // check input GENIE ROOT file
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
    << "\n\n" << "Syntax:" << "\n   test -f input_filename\n";
}
//___________________________________________________________________
