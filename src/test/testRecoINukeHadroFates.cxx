//____________________________________________________________________________
/*!

\program gtestRecoINukeHadroFates

\brief   A simple program to test reconstruction of INTRANUKE/hA hadron fates.

\created March 01, 2009

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <iomanip>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Algorithm/AlgConfigPool.h"
#include "Registry/Registry.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Target.h"
#include "Utils/NuclearUtils.h"
#include "PDG/PDGCodes.h"

using std::string;
using namespace genie;

bool          testCommandLineArgs (int argc);
string        getRootFilename     (int argc, char ** argv);
int           getNEvents          (int argc, char ** argv);
bool          checkRootFilename   (string filename);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- create output tree & file

  int    evt_num;
  int    hadron_pdgcode;
  double hadron_energy;
  int    hadron_fate;

  TFile * output_file = new TFile("hadron_summarytree.root","RECREATE");
  TTree * output_tree = new TTree("output_tree","summary info for events generated with ghAevgen");

  output_tree->Branch("event_num",     &evt_num,        "event_num/I");
  output_tree->Branch("hadron_pdg",    &hadron_pdgcode, "hadron_pdg/I");
  output_tree->Branch("hadron_energy", &hadron_energy,  "hadron_energy/D");
  output_tree->Branch("hadron_fate",   &hadron_fate,    "hadron_fate/I");

  //-- scan the command line arguments and get the input event file
  //-- the number of events to loop over
  string filename = getRootFilename(argc,argv);
  if ( !testCommandLineArgs(argc)   ) return 1;
  if ( !checkRootFilename(filename) ) return 2;

  //-- open the ROOT file and get the TTree & its header
  TFile input_file(filename.c_str(),"READ");
  TTree * input_tree = dynamic_cast <TTree *> (input_file.Get("gtree"));
  NtpMCEventRecord * mcrec = 0;
  input_tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over event tree

  int nevents = getNEvents(argc,argv);
  if(input_tree->GetEntries() < nevents){ nevents = input_tree->GetEntries();}
  LOG("Main", pINFO)
     << "Processing " << nevents << " events";

  for(int i = 0; i < nevents; i++) {
    LOG("Main", pINFO) 
      << "\n\n\n\n--------Determining hadron fates for event "<< i << " -----------";

    input_tree->GetEntry(i);

    EventRecord &  event = *(mcrec->event);
    LOG("Main", pINFO) << event;

    // Get the incoming hadron 
    GHepParticle * incoming_hadron = event.Particle(0);

    // Fill in summary ntuple variables
    evt_num        = i;
    hadron_pdgcode = incoming_hadron->Pdg();
    hadron_energy  = incoming_hadron->P4()->Energy();
    hadron_fate    = (int) utils::intranuke::ReconstructHadronFateHA(&event,0,true); 

    LOG("Main", pNOTICE) 
    << "Event: " << i << ", Hadron code = " << hadron_pdgcode 
    << ", energy = " << hadron_energy << ", fate = " << hadron_fate;

    output_tree->Fill();
    mcrec->Clear();

  }//event loop

  input_file.Close();

  output_file->Write();
  output_file->Close();

  LOG("Main", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
bool testCommandLineArgs(int argc)
{
  if(argc!=5) {
   LOG("Main", pERROR) << "Not enough command line arguments";
   LOG("Main", pINFO)  << "Syntax: gtestGetHadronFates -f root_filename -n nevents";
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
int getNEvents(int argc, char ** argv)
{
  //-- Scan for the number of events to loop over from the command line argument (following -n)
  int nevents=1;
  for(int iarg = 0; iarg < argc-1; iarg++) {
   string argument(argv[iarg]);
   if( argument.compare("-n") == 0 ) nevents = std::atoi(argv[++iarg]);
  }
  return nevents;
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


