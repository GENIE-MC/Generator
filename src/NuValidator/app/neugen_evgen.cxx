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
            ngevgen [-n nev] [-e energy] [-p nupdg] [-t tgtpdg] [-o opt]

         Options :
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

#include <TFile.h>
#include <TTree.h>

#include "Facades/NeuGenWrapper.h"
#include "EventGeneration/EventRecord.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNevents    = 300;        // default n-events to generate
int    gOptPrint      = 3;          // default output (printing) option
double gOptNuEnergy   = 3.0;        // default neutrino energy
int    gOptNuPdgCode  = kPdgNuMu;   // default neutrino PDG code
int    gOptTgtPdgCode = 1056026000; // default target PDG code

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

  //-- initialize an Ntuple Writer
  NtpWriter ntpw;
  ntpw.InitTree("./neugen_events.root");

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
  char * argument = new char[128];

  while( argc>1 && (argv[1][0] == '-'))
  {
    // parse command line argument for number of events
    if (argv[1][1] == 'n') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptNevents = atoi(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptNevents = atoi(argument);
      }
    }
    // parse command line argument for neutrino energy
    if (argv[1][1] == 'e') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptNuEnergy = atof(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptNuEnergy = atof(argument);
      }
    }
    // parse command line argument for neutrino PDG code
    if (argv[1][1] == 'p') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptNuPdgCode = atoi(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptNuPdgCode = atoi(argument);
      }
    }
    // parse command line argument for target pdg code
    if (argv[1][1] == 't') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptTgtPdgCode = atoi(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptTgtPdgCode = atoi(argument);
      }
    }
    // parse command line argument for output (print) options
    if (argv[1][1] == 'o') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptPrint = atoi(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptPrint = atoi(argument);
      }
    }
    argc--;
    argv++;
  }
  delete [] argument;
}
//____________________________________________________________________________

