//____________________________________________________________________________
/*!

\program gvld_sample_scan

\brief   A simple validation program to read-in a GHEP event tree and test 
         whether the generated events obey basic conservation laws

\syntax  gvld_sample_scan 
             -f ghep_event_file 
            [-o output_error_log_file]
            [-n nev1[,nev2]]
            [--check-energy-momentum-conservation]
            [--check-charge-conservation]
            [--check-for-num-of-final-state-nucleons-inconsistent-with-target]
            [--check-for-pseudoparticles-in-final-state]
            [--check-for-off-the-mass-shell-particles-in-final-state]

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
#include <sstream>
#include <fstream>

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
#include "PDG/PDGUtils.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"

using std::ostringstream;
using std::string;
using std::vector;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::endl;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
bool CheckRootFilename  (string filename);

// checks
void CheckEnergyMomentumConservation (void);
void CheckChargeConservation (void);
void CheckForPseudoParticlesInFinState (void);
void CheckForOffMassShellParticlesInFinState (void);

// options
string   gOptInpFilename = "";
string   gOptOutFilename = "";
Long64_t gOptNEvtL = -1;
Long64_t gOptNEvtH = -1;
bool     gOptCheckEnergyMomentumConservation = false;
bool     gOptCheckChargeConservation = false;
bool     gOptCheckForPseudoParticlesInFinState = false;
bool     gOptCheckForOffMassShellParticlesInFinState = false;

Long64_t gFirstEventNum = -1;
Long64_t gLastEventNum  = -1;

TTree *            gEventTree = 0;
NtpMCEventRecord * gMCRec = 0;
ofstream           gErrLog;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  TFile file(gOptInpFilename.c_str(),"READ");

  NtpMCTreeHeader * thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  LOG("gvldtest", pINFO) << "Input tree header: " << *thdr;
  NtpMCFormat_t format = thdr->format;
  if(format != kNFGHEP) {
      LOG("gvldtest", pERROR) 
        << "*** Unsupported event-tree format : "
        << NtpMCFormat::AsString(format);
      file.Close();
      return 3;
  }

  gEventTree = dynamic_cast <TTree *> (file.Get("gtree"));
  gEventTree->SetBranchAddress("gmcrec", &gMCRec);

  Long64_t nev = gEventTree->GetEntries();
  if(gOptNEvtL == -1 && gOptNEvtH == -1) {
    // read all events
    gFirstEventNum = 0;
    gLastEventNum  = nev-1;
  }
  else {
    // read a range of events
    gFirstEventNum = TMath::Max((Long64_t)0,  gOptNEvtL);
    gLastEventNum = TMath::Min(nev-1,        gOptNEvtH);
    if(gLastEventNum - gFirstEventNum < 0) {
      LOG("gevdump", pFATAL) << "Invalid event range";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  }

  
  if(gOptOutFilename.size() == 0) {
     ostringstream logfile;
     logfile << gOptInpFilename << ".errlog";
     gOptOutFilename = logfile.str();
  }
  if(gOptOutFilename != "none") {
     gErrLog.open(gOptOutFilename.c_str());
     gErrLog << "# ..................................................................................." << endl;
     gErrLog << "# Error log for event file " << gOptInpFilename << endl;
     gErrLog << "# ..................................................................................." << endl;
     gErrLog << "# " << endl;
  }

  if (gOptCheckEnergyMomentumConservation) {
	  CheckEnergyMomentumConservation();
  }
  if (gOptCheckChargeConservation) {
          CheckChargeConservation();
  }
  if (gOptCheckForPseudoParticlesInFinState) {
          CheckForPseudoParticlesInFinState();
  }
  if (gOptCheckForOffMassShellParticlesInFinState) {
          CheckForOffMassShellParticlesInFinState();
  }

  if(gOptOutFilename != "none") {
     gErrLog.close();
  }

  return 0;
}
//____________________________________________________________________________
void CheckEnergyMomentumConservation (void)
{
  LOG("gvldtest", pNOTICE) << "Checking energy/momentum conservation...";

  if(gErrLog.is_open()) {
    gErrLog << "# Events failing the energy-momentum conservation test:" << endl;
    gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gvldtest", pINFO) << "Checking event.... " << i;

    double E_init  = 0, E_fin  = 0; // E
    double px_init = 0, px_fin = 0; // px
    double py_init = 0, py_fin = 0; // py
    double pz_init = 0, pz_fin = 0; // pz

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist  = p->Status();

      if(ist == kIStInitialState) 
      {
         E_init  += p->E();
         px_init += p->Px();
         py_init += p->Py();
         pz_init += p->Pz();
       }
       if(ist == kIStStableFinalState || 
          ist == kIStFinalStateNuclearRemnant) 
       {
         E_fin   += p->E();
         px_fin  += p->Px();
         py_fin  += p->Py();
         pz_fin  += p->Pz();
       }
    }//p

    double epsilon = 1E-3; 

    bool E_conserved  = TMath::Abs(E_init  - E_fin)  < epsilon;
    bool px_conserved = TMath::Abs(px_init - px_fin) < epsilon;
    bool py_conserved = TMath::Abs(py_init - py_fin) < epsilon;
    bool pz_conserved = TMath::Abs(pz_init - pz_fin) < epsilon;

    bool ok = E_conserved  && 
              px_conserved &&
              py_conserved &&
              pz_conserved;

    if(!ok) {
       LOG("gvldtest", pERROR) 
         << " ** Energy-momentum non-conservation in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
           gErrLog << i << endl;    
       }
       nerr++;
    }

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
         gErrLog << "none" << endl;    
     }
  }

  LOG("gvldtest", pNOTICE) 
     << "Found " << nerr 
     << " events failing the energy/momentum conservation test";
}
//____________________________________________________________________________
void CheckChargeConservation(void)
{
  LOG("gvldtest", pNOTICE) << "Checking charge conservation...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events failing the charge conservation test:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gvldtest", pINFO) << "Checking event.... " << i;

    // Can't run the test for neutrinos scattered off nuclear targets
    // because of intranuclear rescattering effects and the presence, in the event
    // record, of a charged nuclear remnant pseudo-particle whose charge is not stored.
    // To check charge conservation in the primary interaction, use a sample generated
    // for a free nucleon targets.
    GHepParticle * nucltgt = event.TargetNucleus();
    if (nucltgt) {
      LOG("gvldtest", pINFO)
           << "Event in nuclear target - Skipping test...";
      continue;
    }

    double Q_init  = 0;
    double Q_fin   = 0; 

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist  = p->Status();

      if(ist == kIStInitialState) 
      {
         Q_init  += p->Charge();
       }
       if(ist == kIStStableFinalState)
       {
         Q_fin  += p->Charge();
       }
    }//p

    double epsilon = 1E-3; 
    bool ok = TMath::Abs(Q_init - Q_fin)  < epsilon;
    if(!ok) {
       LOG("gvldtest", pERROR) 
         << " ** Charge non-conservation in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
          gErrLog << i << endl;    
       }
       nerr++;
    }

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gvldtest", pNOTICE) 
     << "Found " << nerr 
     << " events failing the charge conservation test";
}
//____________________________________________________________________________
void CheckForPseudoParticlesInFinState(void)
{
  LOG("gvldtest", pNOTICE) 
      << "Checking for pseudo-particles appearing in final state...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events with pseudo-particles in final state:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gvldtest", pINFO) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    bool ok = true;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist = p->Status();
      if(ist != kIStStableFinalState) continue;
      int pdgc = p->Pdg();      
      if(pdg::IsPseudoParticle(pdgc))
      {
        ok = false;
        break;
      }      
    }//p

    if(!ok) {
       LOG("gvldtest", pERROR) 
         << " ** Pseudo-particle final state particle in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
          gErrLog << i << endl;    
       }
       nerr++;
    }

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gvldtest", pNOTICE) 
     << "Found " << nerr 
     << " events with pseudo-particles in  final state";
}
//____________________________________________________________________________
void CheckForOffMassShellParticlesInFinState(void)
{
  LOG("gvldtest", pNOTICE) 
      << "Checking for off-mass-shell particles appearing in the final state...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events with off-mass-shell particles in final state:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gvldtest", pINFO) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    bool ok = true;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist = p->Status();
      if(ist != kIStStableFinalState) continue;
      if(p->IsOffMassShell())
      {
        ok = false;
        break;
      }      
    }//p

    if(!ok) {
       LOG("gvldtest", pERROR) 
         << " ** Off-mass-shell final state particle in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
          gErrLog << i << endl;    
       }
       nerr++;
    }

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gvldtest", pNOTICE) 
     << "Found " << nerr 
     << " events with off-mass-shell particles in final state";
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);
  
  // get input GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("gvldtest", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gvldtest", pFATAL) 
        << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // get output error log
  if( parser.OptionExists('o') ) {
    LOG("gvldtest", pINFO) << "Reading err log file name";
    gOptOutFilename = parser.ArgAsString('o');
  } 

  // number of events
  if ( parser.OptionExists('n') ) {
    LOG("gevdump", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("gevdump", pFATAL) << "Invalid syntax";
         PrintSyntax();
         gAbortingInErr = true;
         exit(1);
      }
      // read a range of events
      gOptNEvtL = vecn[0];
      gOptNEvtH = vecn[1];
    } else {
      // read single event
      gOptNEvtL = parser.ArgAsLong('n');
      gOptNEvtH = gOptNEvtL;
    }
  } else {
    LOG("gevdump", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvtL = -1;
    gOptNEvtH = -1;
  }

  // checks
  gOptCheckEnergyMomentumConservation =
     parser.OptionExists("check-energy-momentum-conservation");
  gOptCheckChargeConservation =
     parser.OptionExists("check-charge-conservation");
  gOptCheckForPseudoParticlesInFinState = 
     parser.OptionExists("check-for-pseudoparticles-in-final-state");
  gOptCheckForOffMassShellParticlesInFinState = 
     parser.OptionExists("check-for-off-the-mass-shell-particles-in-final-state");
  
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_sample_scan -f sample.root [-n n1[,n2]] [-o errlog] [check names]\n";
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

