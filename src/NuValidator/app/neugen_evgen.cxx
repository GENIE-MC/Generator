//____________________________________________________________________________
/*!

\program ngevgen

\brief   Uses the NeuGEN wrapper to generate neutrino interactions, translates
         the NeuGEN STDHEP common to GENIE's EventRecord and saves the NeuGEN
         events into the GENIE standard output ntuple.

         This program is used for NeuGEN/GENIE cross-comparisons of generated
         event samples using a common ntuple format (GENIE's output ntuple)

         The equivalent GENIE program can be found in src/test/testEvGen.cxx

         Syntax :
            ngevgen [-n nev] -e energy -p nupdg -t tgtpdg [-o opt]

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate
           -e specifies the neutrino energy
           -p specifies the neutrino PDG code
           -t specifies the target PDG code (std format: 1aaazzz000)
           -o specifies the event printing option :
                0 inhibits the event print-out
                1 prints only the original NeuGEN STDHEP common
                2 prints only the translated GENIE EventRecord
                3 prints both NeuGEN and GENIE event records

         Example:
           ngevgen -n 300 -s -e 6.5 -p 14 -t 1056026000 -o 3

           will generate 300 events of muon neutrinos (pdg = 14) on Iron
           (A=56,Z=26) at E = 6.5 GeV and it will be printing both NeuGEN's
           STDHEP common and GENIE's GHEP record

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 20, 2005 at a boring MINOS shift...
*/
//____________________________________________________________________________

#include <sstream>

#include <TFile.h>
#include <TTree.h>

#include "Facades/NeuGenWrapper.h"
#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpWriter.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//Default options (override them using the command line arguments):
int    kDefOptNevents   = 1000; // n-events to generate
int    kDefOptPrint     = 3;    // default output (printing) option
//User-specified options:
int    gOptNevents;             // n-events to generate
int    gOptPrint;               // output (printing) option
double gOptNuEnergy;            // neutrino energy
int    gOptNuPdgCode;           // neutrino PDG code
int    gOptTgtPdgCode;          // target PDG code

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- set custom priority levels for various output streams

  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments
  LOG("ngevgen", pINFO) << "Requested number of events  = " << gOptNevents;
  LOG("ngevgen", pINFO) << "Specified printing option   = " << gOptPrint;
  LOG("ngevgen", pINFO) << "Specified neutrino energy   = " << gOptNuEnergy;
  LOG("ngevgen", pINFO) << "Specified neutrino PDG code = " << gOptNuPdgCode;
  LOG("ngevgen", pINFO) << "Specified target PDG code   = " << gOptTgtPdgCode;

  //-- start NeuGEN
  NeuGenWrapper neugen;

  //-- initialize an Ntuple Writer (build PR rathen that ER ntuples
  //   since event records wrapped from NeuGEN lack the Interaction
  //   summary and it has to be recosntructed)

  NtpMCFormat_t format = kNFPlainRecord;

  ostringstream filename;
  filename << "GNtp" << NtpMCFormat::FilenameTag(format) << "-NeuGEN.root";

  NtpWriter ntpw(format);
  ntpw.InitTree(filename.str());

  //-- initialize NeuGEN event generation
  char * flag = "NOFF";
  int    ctrl = 1;
  neugen.GenControl(flag,ctrl);

  //-- get nuclear target A, Z from input PDG code
  int Z = pdg::IonPdgCodeToZ(gOptTgtPdgCode);
  int A = pdg::IonPdgCodeToA(gOptTgtPdgCode);

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     // generate NeuGEN event and translate it into GENIE's event record
     EventRecord * ev = neugen.GenerateEvent(gOptNuPdgCode,gOptNuEnergy,A,Z);

     if(gOptPrint == 1 || gOptPrint == 3) {
         LOG("ngevgen", pINFO) << "NeuGEN's STDHEP common:";
         neugen.PrintEvent();
         LOG("ngevgen", pINFO) << "\n\n";
     }
     if(gOptPrint == 2 || gOptPrint == 3) {
         LOG("ngevgen", pINFO) << "GENIE's Event Record:";
         LOG("ngevgen", pINFO) << *ev;
     }

     // add event to ntuple
     ntpw.AddEventRecord(ievent++, ev);

     delete ev;
  }
  //-- save the ntuple
  ntpw.SaveTree();

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("ngevgen", pNOTICE) << "Parsing commad line arguments";

  //-- Optional arguments

  //number of events:
  try {
    LOG("ngevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  }
  catch (exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ngevgen", pNOTICE)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }
  //print-out options
  try {
    LOG("ngevgen", pINFO) << "Reading requested print-out option";
    gOptPrint = genie::utils::clap::CmdLineArgAsInt(argc,argv,'o');
  }
  catch (exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ngevgen", pNOTICE)
                          << "Unspecified print option - Using default";
      gOptPrint = kDefOptPrint;
    }
  }

  //-- Required arguments

  //neutrino energy:
  try {
    LOG("ngevgen", pINFO) << "Reading neutrino energy";
    gOptNuEnergy = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'e');
  }
  catch (exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ngevgen", pFATAL) << "Unspecified neutrino energy - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
  //neutrino PDG code:
  try {
    LOG("ngevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  }
  catch (exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ngevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
  //target PDG code:
  try {
    LOG("ngevgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  }
  catch (exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ngevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("ngevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
     << "  ngevgen [-n nev] -e energy -p nupdg -t tgtpdg [-o print_opt]\n";
}
//____________________________________________________________________________

