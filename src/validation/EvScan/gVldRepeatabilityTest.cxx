//____________________________________________________________________________
/*!

\program gvld_repeatability_test

\brief   Tests the repeatability of GENIE MC sample generation.
         The program compares two input event files which must be generated,
         in separate processes, with the exact same inputs and with the same 
         random number seed number.

         Syntax:
           shell$ gvld_repeatability_test
                    --first-sample /path/to/file.root
                    --second-sample /path/to/file.root
                   [--event-record-print-level level]
                   [--add-event-printout-in-error-log]
                   [--max-num-of-errors-shown n]
                   [-o output_error_log_file]
                   [-n n1[,n2]] 

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created April 23, 2013

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <fstream>

#include <TFile.h>
#include <TTree.h>

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/RunOpt.h"

using std::string;
using std::endl;
using std::ostream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
void GetEventRange      (Long64_t nev, Long64_t & n1, Long64_t & n2);

Long64_t gOptNEvtL;
Long64_t gOptNEvtH;
int      gOptMaxNumErrs = -1; 
bool     gOptAddEventPrintoutInErrLog = false;
string   gOptInpFilename1;
string   gOptInpFilename2;
string   gOptOutFilename = "";

const double epsilon = 1E-12;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // set print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  //
  // open the ROOT files and get the GHEP event trees
  //

  TFile file1(gOptInpFilename1.c_str(),"READ");
  TTree * ghep_tree1 = 
     dynamic_cast <TTree *> (file1.Get("gtree"));
  if(!ghep_tree1) {
    LOG("gvldtest", pFATAL) 
        << "No GHEP event tree in input file: " << gOptInpFilename1;
    gAbortingInErr=true;
    exit(1);
  }
  Long64_t nev1 = ghep_tree1->GetEntries();
  LOG("gvldtest", pFATAL) 
     << "The first input GHEP event tree has " << nev1 
     << ((nev1==1) ? " entry." : " entries.");

  TFile file2(gOptInpFilename2.c_str(),"READ");
  TTree * ghep_tree2 = 
     dynamic_cast <TTree *> (file2.Get("gtree"));
  if(!ghep_tree2) {
    LOG("gvldtest", pFATAL) 
        << "No GHEP event tree in input file: " << gOptInpFilename2;
    gAbortingInErr=true;
    exit(1);
  }
  Long64_t nev2 = ghep_tree2->GetEntries();
  LOG("gvldtest", pFATAL) 
     << "The second input GHEP event tree has " << nev2 
     << ((nev2==1) ? " entry." : " entries.");

  assert(nev1==nev2);

  //  
  // set branch addresses
  //

  NtpMCEventRecord * mcrec1 = 0;
  ghep_tree1->SetBranchAddress("gmcrec", &mcrec1);
  NtpMCEventRecord * mcrec2 = 0;
  ghep_tree2->SetBranchAddress("gmcrec", &mcrec2);

  //
  // output error log
  //

  ofstream errLog;
  if(gOptOutFilename.size() > 0) {
     errLog.open(gOptOutFilename.c_str());
     errLog << "# ..................................................................................." << endl;
     errLog << "# Comparison of " << endl;
     errLog << "#  - " << gOptInpFilename1 << endl;
     errLog << "#  - " << gOptInpFilename2 << endl;
     errLog << "# ..................................................................................." << endl;
     errLog << "# " << endl;
  }

  //
  // event loop
  //

  int nerr = 0;

  Long64_t n1,n2;
  GetEventRange(nev1,n1,n2);
  for(Long64_t i = n1; i <= n2; i++) {

    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    // get events
    ghep_tree1->GetEntry(i);
    EventRecord & event1 = *(mcrec1->event);
    ghep_tree2->GetEntry(i);
    EventRecord & event2 = *(mcrec2->event);

    //
    // first, check generated info with event wide-scope
    //
    bool same = 
      ( TMath::Abs(event1.Weight()      - event2.Weight())      < epsilon ) &&
      ( TMath::Abs(event1.Probability() - event2.Probability()) < epsilon ) &&
      ( TMath::Abs(event1.XSec()        - event2.XSec())        < epsilon ) &&
      ( TMath::Abs(event1.DiffXSec()    - event2.DiffXSec())    < epsilon ) &&
      ( TMath::Abs(event1.Vertex()->X() - event2.Vertex()->X()) < epsilon ) &&
      ( TMath::Abs(event1.Vertex()->Y() - event2.Vertex()->Y()) < epsilon ) &&
      ( TMath::Abs(event1.Vertex()->Z() - event2.Vertex()->Z()) < epsilon ) &&
      ( TMath::Abs(event1.Vertex()->T() - event2.Vertex()->T()) < epsilon );

    if(!same) {
       nerr++;
       LOG("gvldtest", pERROR)
            << "Difference seen in event: " << i
            << "\n"
            << event1
            << event2;
       if(errLog.is_open()) {
           errLog << "Difference seen in event: " << i << endl;
           if(gOptAddEventPrintoutInErrLog) {
               errLog << event1;
               errLog << event2;
           }
       }
       mcrec1->Clear();
       mcrec2->Clear();
       continue;
    }

    //
    // then, compare all initial / intermediate and final state particles one by one.
    //
    Int_t np1 = event1.GetEntries();
    Int_t np2 = event2.GetEntries();
    if(np1 != np2) {
    }
    bool same_particle_list = true;
    for(int j = 0; j < np1; j++) {
       GHepParticle * p1 = event1.Particle(j);
       GHepParticle * p2 = event2.Particle(j);
       bool same_particle =
         ( p1->Pdg()           == p2->Pdg()           ) &&
         ( p1->Status()        == p2->Status()        ) &&
         ( p1->RescatterCode() == p2->RescatterCode() ) &&
         ( p1->FirstMother()   == p2->FirstMother()   ) &&
         ( p1->LastMother()    == p2->LastMother()    ) &&
         ( p1->FirstDaughter() == p2->FirstDaughter() ) &&
         ( p1->LastDaughter()  == p2->LastDaughter()  ) &&
         ( TMath::Abs(p1->Px() - p2->Px()) < epsilon  ) &&
         ( TMath::Abs(p1->Py() - p2->Py()) < epsilon  ) &&
         ( TMath::Abs(p1->Pz() - p2->Pz()) < epsilon  ) &&
         ( TMath::Abs(p1->E()  - p2->E())  < epsilon  ) &&
         ( TMath::Abs(p1->Vx() - p2->Vx()) < epsilon  ) &&
         ( TMath::Abs(p1->Vy() - p2->Vy()) < epsilon  ) &&
         ( TMath::Abs(p1->Vz() - p2->Vz()) < epsilon  ) &&
         ( TMath::Abs(p1->Vt() - p2->Vt()) < epsilon  ) &&
         ( TMath::Abs(p1->PolzPolarAngle()   - p2->PolzPolarAngle()  ) < epsilon  ) &&
         ( TMath::Abs(p1->PolzAzimuthAngle() - p2->PolzAzimuthAngle()) < epsilon  );
      if(!same_particle) {
         same_particle_list = false;
         break;
      }
    }//j
    if(!same_particle_list) {
       nerr++;
       LOG("gvldtest", pERROR)
            << "Difference seen in event: " << i
            << "\n"
            << event1
            << event2;
       if(errLog.is_open()) {
           errLog << "Difference seen in event: " << i << endl;
           if(gOptAddEventPrintoutInErrLog) {
               errLog << event1;
               errLog << event2;
           }
       }
       mcrec1->Clear();
       mcrec2->Clear();
       continue;
    }

    mcrec1->Clear();
    mcrec2->Clear();

  }//i

  // clean-up
  file1.Close();
  file2.Close();

  if(nerr == 0) {
    LOG("gvldtest", pNOTICE) << "The files are identical";
    if(errLog.is_open()) {
       errLog << "The files are identical" << endl;
    }
  }

  LOG("gevdump", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetEventRange(Long64_t nev, Long64_t & n1, Long64_t & n2)
{
  if(gOptNEvtL == -1 && gOptNEvtH == -1) {
    // read all events
    n1=0;
    n2=nev-1;
  }
  else {
    // read a range of events
    n1 = TMath::Max((Long64_t)0,  gOptNEvtL);
    n2 = TMath::Min(nev-1,        gOptNEvtH);
    if(n2-n1 <0) {
      LOG("gevdump", pFATAL) << "Invalid event range";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    } 
  }
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevdump", pINFO) << "*** Parsing command line arguments";

  // Common run options.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // get GENIE event samples
  if ( parser.OptionExists("first-sample") ) {
    LOG("gevdump", pINFO) << "Reading 1st event sample filename";
    gOptInpFilename1 = parser.ArgAsString("first-sample");
  } else {
    LOG("gevdump", pFATAL) 
       << "Unspecified 1st sample - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }
  if ( parser.OptionExists("second-sample") ) {
    LOG("gevdump", pINFO) << "Reading 2nd event sample filename";
    gOptInpFilename2 = parser.ArgAsString("second-sample");
  } else {
    LOG("gevdump", pFATAL) 
       << "Unspecified 2nd sample - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
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

  // get output error log
  if( parser.OptionExists('o') ) {
    LOG("gvldtest", pINFO) << "Reading err log file name";
    gOptOutFilename = parser.ArgAsString('o');
  } 


  gOptAddEventPrintoutInErrLog =
     parser.OptionExists("add-event-printout-in-error-log");

  if(parser.OptionExists("max-num-of-errors-shown")) {
     gOptMaxNumErrs = parser.ArgAsInt("max-num-of-errors-shown");
     gOptMaxNumErrs = TMath::Max(1,gOptMaxNumErrs);
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevdump", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevdump -f sample.root [-n n1[,n2]] [--event-record-print-level]\n";
}
//_________________________________________________________________________________
