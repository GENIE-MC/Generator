//____________________________________________________________________________
/*!

\program gvld_particle_test

\brief   

\syntax  gvld_particle_test -f ghep_event_tree_filename

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 13, 2008

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <map>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGLibrary.h"
#include "Messenger/Messenger.h"

using std::string;
using namespace genie;

bool   testCommandLineArgs(int argc);
string getRootFilename    (int argc, char ** argv);
bool   checkRootFilename  (string filename);

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //
  // ** initialization
  //

  // scan the command line arguments and get the ROOT filename
  string filename = getRootFilename(argc,argv);

  if ( !testCommandLineArgs(argc)   ) return 1;
  if ( !checkRootFilename(filename) ) return 2;

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(filename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("gvldtest", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  if(format != kNFGHEP) {
      LOG("gvldtest", pERROR) 
        << "*** Unsupported event-tree format : "
        << NtpMCFormat::AsString(format);
      file.Close();
      return 3;
  }

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //
  // ** event loop / analysis
  //

  map<int, int> fs_particles;

  //-- loop over events & test conservation laws
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gvldtest", pNOTICE) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     if(p->Status() == kIStStableFinalState) {
       int pdgc = p->Pdg();
       fs_particles[pdgc]++;
     }
    
    } // particle loop

    mcrec->Clear();

  } // event loop

  file.Close();


  //
  // ** Reporting
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();



  LOG("gvldtest", pINFO)  << "Done!";

  return 0;
}
//____________________________________________________________________________
bool testCommandLineArgs(int argc)
{
  if(argc!=3) {
   LOG("gvldtest", pERROR) 
       << "\n"
       << "*** Not enough command line arguments"
       << "    Syntax: gvld_conservation_test -f ghep_event_tree_filename";
   return false;
  }
  return true;
}
//____________________________________________________________________________
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
//____________________________________________________________________________
bool checkRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gvldtest", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//____________________________________________________________________________

