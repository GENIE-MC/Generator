//____________________________________________________________________________
/*!

\program gvld_particle_test

\brief   

\syntax  gvld_particle_test -f ghep_event_tree_filename [-n nevents]

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
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

string gOptInpFileName;
int    gOptNEvt;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);

  //
  // ** initialization
  //

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFileName.c_str(),"READ");

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

  map<int, int> all_particles;   // all
  map<int, int> fs_particles;    // final state particles
  map<int, int> primh_particles; // primary hadronic system particles
  map<int, int> dec_particles;   // particles marked as decayed

  long int all_n   = 0;
  long int fs_n    = 0;
  long int primh_n = 0;
  long int dec_n   = 0;

  //-- event loop

  int nev = (gOptNEvt > 0) ? 
                 TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
                 (int) tree->GetEntries();

  for(int i = 0; i< nev; i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gvldtest", pNOTICE) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     int pdgc = p->Pdg();

     all_particles[pdgc]++;
     all_n++;
    
     if (p->Status() == kIStStableFinalState) 
     {
       fs_particles[pdgc]++;
       fs_n++;
     }
     else if (p->Status() == kIStDecayedState) 
     {
       dec_particles[pdgc]++;
       dec_n++;
     }
     else if (p->Status() == kIStHadronInTheNucleus) 
     {
       primh_particles[pdgc]++;
       primh_n++;
     }    

    } // particle loop

    mcrec->Clear();

  } // event loop

  file.Close();


  //
  // ** Reporting
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();

  map<int, int>::const_iterator it;; 

  LOG("gvldtest", pINFO)  << "*** All particles:";
  for(it = all_particles.begin(); it != all_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (all_n>0) ? (double)n / (double)all_n : 0.;
     LOG("gvldtest", pINFO)  
       << " -> " << pdglib->Find(pdgc)->GetName() << " : " 
       << n << " (" << 100*r << " %)";
  }
  LOG("gvldtest", pINFO)  << "*** Final state particles:";
  for(it = fs_particles.begin(); it != fs_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (fs_n>0) ? (double)n / (double)fs_n : 0.;
     LOG("gvldtest", pINFO)  
       << " -> " << pdglib->Find(pdgc)->GetName() << " : " 
       << n << " (" << 100*r << " %)";
  }
  LOG("gvldtest", pINFO)  << "*** Decayed particles:";
  for(it = dec_particles.begin(); it != dec_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (dec_n>0) ? (double)n / (double)dec_n : 0.;
     LOG("gvldtest", pINFO)  
       << " -> " << pdglib->Find(pdgc)->GetName() << " : " 
       << n << " (" << 100*r << " %)";
  }
  LOG("gvldtest", pINFO)  << "*** Primary hadronic system:";
  for(it = primh_particles.begin(); it != primh_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (primh_n > 0) ? (double)n / (double)primh_n : 0;
     LOG("gvldtest", pINFO)  
       << " -> " << pdglib->Find(pdgc)->GetName() << " : " 
       << n << " (" << 100*r << " %)";
  }

  LOG("gvldtest", pINFO)  << "Done!";

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  //get input filename
  try {
    gOptInpFileName = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gvldtest", pFATAL)
        << "Unspecified input filename - Exiting";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  }

  // number of events:
  try {
    gOptNEvt = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptNEvt = -1;
    }
  }   
}
//____________________________________________________________________________
void PrintSyntax(void)
{
   LOG("gvldtest", pERROR) 
       << "\n"
       << "Syntax:\n"
       << "gvld_particle_test -f ghep_file [-n nevt]";
}
//____________________________________________________________________________
