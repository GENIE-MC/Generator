//____________________________________________________________________________
/*!

\program gvld_conservation_test

\brief   A simple validation program to read-in a GHEP event tree and test 
         whether the generated events obey basic conservation laws

\syntax  gvld_conservation_test -f ghep_event_tree_filename

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 13, 2008

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

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

  double E_init  = 0, E_fin  = 0; // E
  double px_init = 0, px_fin = 0; // px
  double py_init = 0, py_fin = 0; // py
  double pz_init = 0, pz_fin = 0; // pz
  double Q_init  = 0, Q_fin  = 0; // charge
  double s_init  = 0, s_fin  = 0; // strangeness

  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- loop over events & test conservation laws
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gvldtest", pNOTICE) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     if(p->Status() == kIStInitialState) {
       E_init   += p->E();
       px_init  += p->Px();
       py_init  += p->Py();
       pz_init  += p->Pz();
       Q_init   += p->Charge();
       s_init   += pdglib->Find(p->Pdg())->Strangeness();
     }

     if(p->Status() == kIStStableFinalState) {
       E_fin   += p->E();
       px_fin  += p->Px();
       py_fin  += p->Py();
       pz_fin  += p->Pz();
       Q_fin   += p->Charge();
       s_fin   += pdglib->Find(p->Pdg())->Strangeness();
     }
    
    } // particle loop

    double epsilon = 1E-3; 

    bool E_conserved  = TMath::Abs(E_init  - E_fin)  < epsilon;
    bool px_conserved = TMath::Abs(px_init - px_fin) < epsilon;
    bool py_conserved = TMath::Abs(py_init - py_fin) < epsilon;
    bool pz_conserved = TMath::Abs(pz_init - pz_fin) < epsilon;
    bool Q_conserved  = TMath::Abs(Q_init  - Q_fin)  < epsilon;
    bool s_conserved  = TMath::Abs(s_init  - s_fin)  < epsilon;

    bool ok = E_conserved  && 
              px_conserved &&
              py_conserved &&
              pz_conserved &&
              Q_conserved  &&
              s_conserved;

    if(!ok) {
       if(!E_conserved) {
          LOG("gvldtest", pERROR) << "** E is not conserved at this event";
       }
       if(!px_conserved) {
          LOG("gvldtest", pERROR) << "** Px is not conserved at this event";
       }
       if(!py_conserved) {
          LOG("gvldtest", pERROR) << "** Py is not conserved at this event";
       }
       if(!pz_conserved) {
          LOG("gvldtest", pERROR) << "** Pz is not conserved at this event";
       }
       if(!Q_conserved) {
          LOG("gvldtest", pERROR) << "** Q is not conserved at this event";
       }
       if(!s_conserved) {
          LOG("gvldtest", pERROR) << "** s is not conserved at this event";
       }

       LOG("gvldtest", pERROR) << rec_header;
       LOG("gvldtest", pERROR) << event;
    } 

    mcrec->Clear();

  } // event loop

  file.Close();

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

