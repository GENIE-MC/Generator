//____________________________________________________________________________
/*!

\program gtestRwNuXSecHelperLoop

\brief   A simple program to test the GReWeightNuXSecHelper class.
         The program will re-weight events generated/stored at some point
         during the generator evolution to take into account the effect of 
         the uncertainty at a set of physics parameters.
         So reweighting happens for a given set of physics model but between
         a range of input physics parameters. 
	 Most continuous (user-space) physics parameters listed in 
         $GENIE/config/UserPhysicsOptions.xml can be included in the
         parameter loop.

\syntax  gtestRwNuXSecHelperLoop -f filename
         where the filename points to a ROOT file with a GENIE event tree

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 11, 2007

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>

#include "Algorithm/AlgId.h"
#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
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

  // Access the algorithm factory and the configuration pool
  AlgFactory *    algf = AlgFactory::Instance();
  AlgConfigPool * algc = AlgConfigPool::Instance();

  // Get the physics parameters that GENIE developers find meaningfull
  // for a user to tweak
  Registry * user_conf = algc->GlobalParameterList();
  user_conf->UnLock();

  // Loop over physics configuration parameters, tweak their value
  // and reconfigure the entire system of instantiated GENIE algorithms
  // (a simple case here - few values of a single parameter)
  const int npv = 3;
  double MaQEL[npv] = { 1.032, 1.111, 1.211 };

  for(int j = 0; j < npv; j++) 
  {
    user_conf->Set("QEL-Ma", MaQEL[j]);
    algf->ForceReconfiguration();

    LOG("test", pINFO) << "User options / current " << *user_conf;

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

    } // event loop
  } // physics parameter loop

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
    << "\n gtestRwNuXSecHelperLoop -f ghep_file\n";
}
//___________________________________________________________________
