//____________________________________________________________________________
/*!

\program ghAgen

\brief   Generates hadron + nucleus interactions using GENIE's INTRANUKE

         Syntax :
           ghAgen [-n nev] -p hadron_pdg -t tgt_pdg [-r run#] -k KE

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate (default: 10000)
           -p specifies the incoming hadron PDG code 
           -t specifies the nuclear target PDG code (10LZZZAAAI)
           -r specifies the MC run number (default: 0)
           -k specifies the incoming hadron's kinetic energy (in GeV)

         Example:
           ghAgen -n 100000 -p 211 -t 1000260560 -k 0.165
           will generate 100k pi^{+} + Fe56 events at E(pi+)=165 MeV

\author  Minsuk Kim and Steve Dytman
         University of Pittsburgh

\version 1.2

\created May 1, 2007

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TNtuple.h>

#include <iomanip>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

//#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
//#include "Conventions/Controls.h"
//#include "Conventions/Units.h"
//#include "Conventions/Constants.h"
#include "PDG/PDGLibrary.h"
//#include "PDG/PDGCodeList.h"
//#include "Utils/PrintUtils.h"
//#include "Utils/NuclearUtils.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/Intranuke.h"
//#include "HadronTransport/INukeHadroData.h"
//#include "HadronTransport/INukeHadroFates.h"
#include "EVGCore/EventRecordVisitorI.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//using namespace genie::utils;
//using namespace genie::constants;
//using namespace genie::controls;

using std::cout;
using std::endl;
using std::ios;
using std::setw;

//Default options 
int     kDefOptNevents   = 10000;   // n-events to generate
Long_t  kDefOptRunNu     = 0;       // default run number

//User-specified options:
int           gOptNevents;           // n-events to generate
int           gOptTgtPdgCode;        // target PDG code
Long_t        gOptRunNu;             // run number
int           gOptInputPdgCode;      // rescattering particle PDG code as a input
double        gOptInputKE;           // This is KE = E - M. So E = M + KE

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * intranuke = 
            dynamic_cast<const EventRecordVisitorI *> (
                     algf->GetAlgorithm("genie::Intranuke","hA-TestMode"));

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kNFEventRecord, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get the pdg library
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- input 4-momenta & dummy vtx
  double mh  = pdglib->Find(gOptInputPdgCode)->Mass();
  double M   = pdglib->Find(gOptTgtPdgCode)->Mass();

  double Eh  = mh + gOptInputKE;
  double pzh = TMath::Sqrt(TMath::Max(0.,Eh*Eh-mh*mh));

  TLorentzVector p4h   (0.,0.,pzh,Eh);
  TLorentzVector p4tgt (0.,0.,0.,M);
  TLorentzVector x4null(0.,0.,0.,0.);

  //-- generate events / print the GHEP record / add it to the event tree
  int ievent = 0;
  while (ievent < gOptNevents) {
      LOG("ghAgen", pINFO) << " *** Generating event............ " << ievent;
      
      // insert initial state
      EventRecord * evrec = new EventRecord();
      Interaction * interaction = new Interaction;
      evrec->AttachSummary(interaction);
            
      evrec->AddParticle(gOptInputPdgCode,kIStInitialState,-1,-1,-1,-1,p4h,  x4null);
      evrec->AddParticle(gOptTgtPdgCode  ,kIStInitialState,-1,-1,-1,-1,p4tgt,x4null);
      
      // generate h+A eventw
      intranuke->ProcessEventRecord(evrec);
  
      // print generated event    
      LOG("ghAgen", pINFO) << *evrec;
  
      // add event at the output ntuple
      ntpw.AddEventRecord(ievent, evrec);
      
      // refresh the mc job monitor
      mcjmonitor.Update(ievent,evrec);
      
      ievent++;
      delete evrec;
      
  } // end loop events
  
  //-- save the generated MC events
  ntpw.Save();

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("ghAgen", pNOTICE) << "Parsing command line arguments";

  //number of events:
  try {
    LOG("ghAgen", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAgen", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  //run number:
  try {
    LOG("ghAgen", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAgen", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  // incoming hadron PDG code:
  try {
    LOG("ghAgen", pINFO) << "Reading rescattering particle PDG code";
    gOptInputPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
      if(!e.ArgumentFound()) {
      LOG("ghAgen", pFATAL) << "Unspecified PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //target PDG code:
  try {
    LOG("ghAgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // incoming hadron kinetic energy:
  try {
    LOG("ghAgen", pINFO) << "Reading rescattering particle KE energy";
    string ke = genie::utils::clap::CmdLineArgAsString(argc,argv,'k');
    gOptInputKE = atof(ke.c_str());
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAgen", pFATAL) << "Unspecified KE - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  LOG("ghAgen", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("ghAgen", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("ghAgen", pINFO) << "Incoming hadron PDG code   = " << gOptInputPdgCode;
  LOG("ghAgen", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("ghAgen", pINFO) << "Hadron input KE            = " << gOptInputKE;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("ghAgen", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "   ghAgen [-n nev] -p hadron_pdg -t tgt_pdg [-r run] [-a R0] -k KE"
    << "\n\n";
}
//____________________________________________________________________________
