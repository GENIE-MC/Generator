//____________________________________________________________________________
/*!

\program gvld_sample_scan

\brief   A simple validation program to read-in a GHEP event tree and test 
         whether the generated events obey basic conservation laws

\syntax  gvld_samplev_scan -f ghep_event_file [-n nevents]

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 13, 2008

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <vector>
#include <iomanip>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
bool CheckRootFilename  (string filename);
int  Strangeness        (int pdgc);

int    gOptNEvt;
string gOptInpFilename;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

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

  PDGLibrary * pdglib = PDGLibrary::Instance();

  //
  map<int, int> all_particles;   // all
  map<int, int> fs_particles;    // final state particles
  map<int, int> primh_particles; // primary hadronic system particles
  map<int, int> dec_particles;   // particles marked as decayed

  long int all_n   = 0;
  long int fs_n    = 0;
  long int primh_n = 0;
  long int dec_n   = 0;

  double E_init  = 0, E_fin  = 0; // E
  double px_init = 0, px_fin = 0; // px
  double py_init = 0, py_fin = 0; // py
  double pz_init = 0, pz_fin = 0; // pz
  double Q_init  = 0, Q_fin  = 0; // charge
  double s_init  = 0, s_fin  = 0; // strangeness

  vector<int> evtnum;

  int nev = (gOptNEvt > 0) ?
                 TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
                 (int) tree->GetEntries();

  // Event loop

  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gvldtest", pNOTICE) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     int          pdgc = p->Pdg();
     GHepStatus_t ist  = p->Status();

     //
     // Check energy & momentum, charge and strangeness conservation
     //

     if(ist == kIStInitialState) {
       E_init  += p->E();
       px_init += p->Px();
       py_init += p->Py();
       pz_init += p->Pz();
       Q_init  += p->Charge();
     //s_init  += pdglib->Find(p->Pdg())->Strangeness();
       s_init  += Strangeness(pdgc);

     }
     if(ist == kIStStableFinalState || 
        ist == kIStFinalStateNuclearRemnant) {
       E_fin   += p->E();
       px_fin  += p->Px();
       py_fin  += p->Py();
       pz_fin  += p->Pz();
       Q_fin   += p->Charge();
     //s_fin   += pdglib->Find(p->Pdg())->Strangeness();
       s_fin   += Strangeness(pdgc);
     }


     //
     // Check particle frequencies
     //

     all_particles[pdgc]++;
     all_n++;
  
     if (ist == kIStStableFinalState) {
       fs_particles[pdgc]++;
       fs_n++;
     }
     else if (ist == kIStDecayedState) {
       dec_particles[pdgc]++;
       dec_n++;
     }
     else if (ist == kIStHadronInTheNucleus) {
       primh_particles[pdgc]++;
       primh_n++;
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

       evtnum.push_back(i);
    } 

    mcrec->Clear();

  } // event loop

  file.Close();

  //
  // ** Reporting
  //

  LOG("gvldtest", pINFO)  << "Found " << evtnum.size() << " problems";
  for(vector<int>::const_iterator it = evtnum.begin();
      it != evtnum.end(); ++it) {
     LOG("gvldtest", pINFO)  << "Look at event: " << *it;    
  }

 map<int, int>::const_iterator it;

  LOG("gvldtest", pINFO)  << "*** All particles:";
  for(it = all_particles.begin(); it != all_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (all_n>0) ? (double)n / (double)all_n : 0.;
     LOG("gvldtest", pINFO)
       << setfill(' ') << setw(15) 
       << pdglib->Find(pdgc)->GetName() << " :"
       << setfill('.') << setw(15) 
       << n << " (" << 100*r << " %)";
  }
  LOG("gvldtest", pINFO)  << "*** Final state particles:";
  for(it = fs_particles.begin(); it != fs_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (fs_n>0) ? (double)n / (double)fs_n : 0.;
     LOG("gvldtest", pINFO)
       << setfill(' ') << setw(15) 
       << pdglib->Find(pdgc)->GetName() << " :"
       << setfill('.') << setw(15) 
       << n << " (" << 100*r << " %)";
  }
  LOG("gvldtest", pINFO)  << "*** Decayed particles:";
  for(it = dec_particles.begin(); it != dec_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (dec_n>0) ? (double)n / (double)dec_n : 0.;
     LOG("gvldtest", pINFO)
       << setfill(' ') << setw(15) 
       << pdglib->Find(pdgc)->GetName() << " :"
       << setfill('.') << setw(15) 
       << n << " (" << 100*r << " %)";
  } 

  LOG("gvldtest", pINFO)  << "*** Primary hadronic system:";
  for(it = primh_particles.begin(); it != primh_particles.end(); ++it) {
     int pdgc = it->first;
     int n    = it->second;
     double r = (primh_n > 0) ? (double)n / (double)primh_n : 0;
     LOG("gvldtest", pINFO)
       << setfill(' ') << setw(15) 
       << pdglib->Find(pdgc)->GetName() << " :"
       << setfill('.') << setw(15) 
       << n << " (" << 100*r << " %)";
  }

  LOG("gvldtest", pINFO)  << "Done!";

  return 0;
}
//____________________________________________________________________________
int Strangeness(int pdgc)
{
  switch (pdgc) {

   case ( kPdgKP         ) : // K^+
   case ( kPdgK0         ) : // K^0
   case ( kPdgDPs        ) : // Ds^+
   case ( kPdgAntiLambda ) : // \bar{Lambda^0}
   case ( kPdgAntiSigmaP ) : // \bar{Sigma+}
   case ( kPdgAntiSigma0 ) : // \bar{Sigma0}
   case ( kPdgAntiSigmaM ) : // \bar{Sigma-}
   
    return 1;
    break;

   case ( kPdgKM         ) : // K^-
   case ( kPdgAntiK0     ) : // \bar{K^0}
   case ( kPdgDMs        ) : // Ds^-
   case ( kPdgLambda     ) : // Lambda^0
   case ( kPdgSigmaP     ) : // Sigma+
   case ( kPdgSigma0     ) : // Sigma0
   case ( kPdgSigmaM     ) : // Sigma-

    return -1;
    break;

   default:
    return 0;
  }

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // number of events:
  if( parser.OptionExists('n') ) {
    LOG("gvldtest", pINFO) << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("gvldtest", pINFO)
        << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
  
  // get GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("gvldtest", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gvldtest", pFATAL) 
        << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_sample_scan -f sample.root [-n nev] \n";
}
//_________________________________________________________________________________
bool CheckRootFilename(string filename)
{
  if(filename.size() == 0) return false;
    
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gvldtest", pERROR)  
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________

