//____________________________________________________________________________
/*!

\program ngevgen

\brief   Uses the NeuGEN wrapper to generate neutrino interactions, translates
         the NeuGEN STDHEP common to GENIE's EventRecord and saves the NeuGEN
         events into the GENIE standard output ntuple.

         This program is used for NeuGEN / GENIE cross-comparisons of generated
         event samples using a common ntuple format (GENIE's output ntuple)

         The equivalent GENIE program can be found in %GENIE/src/test/testEvGen

         Syntax : neugen_evgen [-n number-of-events] [-p opt]

         -n specifies the number of events to generate 
         -p specifies the event printing option :
                  0 inhibits the event print-out
                  1 prints only the original NeuGEN STDHEP common
                  2 prints only the translated GENIE EventRecord
                  3 prints both NeuGEN and GENIE event records

         Defaults: -n 100 -p 3
         Example : ngevgen -n 300 -s
         
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

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

int gOptNevents = 100; // default n-events to generate
int gOptPrint   = 3;   // default printing option

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- set custom priority levels for various output streams

  GetCommandLineArgs(argc,argv);
  
  //-- print the options you got from command line arguments
  LOG("ngevgen", pINFO) << "Requested number of events = " << gOptNevents;
  LOG("ngevgen", pINFO) << "Specified printing option  = " << gOptPrint;
  
  //-- start NeuGEN
  NeuGenWrapper neugen;

  //-- initialize an Ntuple Writer
  NtpWriter ntpw;
  ntpw.InitTree("./neugen_events.root");
  
  //-- initialize NeuGEN event generation
  char * flag = "NOFF";
  int    ctrl = 1;  
  neugen.GenControl(flag,ctrl);
  
  //-- generate events / print the GHEP record / add it to the ntuple  
  int ievent = 0;
  while ( ievent < gOptNevents) {

     // generate NeuGEN event and translate it into GENIE's event record
     EventRecord * ev_rec = neugen.GenerateEvent(kPdgNuMu,6.782,56,26);

     if(gOptPrint == 1 || gOptPrint == 3) {
         LOG("ngevgen", pINFO) << "NeuGEN's STDHEP common:";
         neugen.PrintEvent();
         LOG("ngevgen", pINFO) << "\n\n";
     }     
     if(gOptPrint == 2 || gOptPrint == 3) {
         LOG("ngevgen", pINFO) << "GENIE's Event Record:";
         LOG("ngevgen", pINFO) << *ev_rec;
     }
     
     // add event to ntuple
     ntpw.AddEventRecord(ievent++, ev_rec);

     delete ev_rec;
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
    if (argv[1][1] == 'p') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptNevents = atoi(argument);
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

