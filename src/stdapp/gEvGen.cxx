//____________________________________________________________________________
/*!

\program gEvGen

\brief   Example program driving GENIE event generation modules

         Syntax :
           gEvGen [-n nev] [-s] [-e energy] [-p nupdg] [-t tgtpdg] [-f format]

         Options :
           -n specifies the number of events to generate
           -s turns on cross section spline building at job initialization
           -e specifies the neutrino energy
           -p specifies the neutrino PDG code
           -t specifies the target PDG code (std format: 1aaazzz000)
           -f specifies the output TTree format. If set to 0 will create a
              single-branch TTree with NtpMCPlainRecord objects in its leaves,
              while if set to 1 it will have NtpMCEventRecord objects in
              its leaves (see the Ntuple package for descriptions of the ntuple
              records and their intended usage). Default options is 1.

         Example:
           gEvGev -n 300 -s -e 6.5 -p 14 -t 1056026000

           will build cross section splines during the initialization step,
           and will generate 300 events of muon neutrinos (pdg = 14) on Iron
           (A=56,Z=26) at E = 6.5 GeV.

         Other control options:

         You can further control the program behaviour by setting the GEVGL,
         GSPLOAD and GSPSAVE environmental variables.

           - Set the GEVGL environmental variable to contol the list of event
             generator objects that get loaded (at job initialization, the
             EventGeneratorListAssember assembles the list of all requested
             event generator objects -meaning: all concrete implementations of
             the EventGeneratorI interface that know how to generate various
             classes of events-). You can see the possible GEVGL options by
             looking up the names of the EventGeneratorListAssembler configuration
             sets at $GENIE/config/event_generator_list_assembler.xml
             If GEVGL is not set, the 'Default' generator list is used.

           - If you specify the -s option you can set the GSPSAVE env.variable
             to specify an XML filename for saving the computed xsec splines.
             If GSPSAVE is not set the computed splines will not be saved.

           - If you specify the -s option you can set the GSPLOAD env.variable
             to specify an XML filename from which to load xsec splines at
             initialization. GENIE will load the splines and then go on and
             compute only the xsec splines that it needs but were not loaded.
             If GSPLOAD is not set, GENIE will not attempt to load any splines
             and will compute them all.

           - You can set the GMSGCONF env.variable to point to a messenger
             XML configuration file (following the syntax of the default one
             that can be found in $GENIE/config/messenger.xml) and modify the
             verbosity of GENIE output. Both the default and your messenger
             configuration will be read but yours will take precedence in
             case of clashing priority for the same stream.

             Examples:
             By setting the following env.vars you ask GENIE to generate QEL
             events only, and load the cross section splines from splines.xml
             and also to add the priority levels in mymsg.xml on top of the
             default ones.
             shell% export GEVGL=QEL
             shell% export GSPLOAD=/home/me/GENIE/mydata/splines.xml
             shell% export GMSGCONF=/home/me/GENIE/mydata/mymsg.xml

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 05, 2004
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "Conventions/XmlParserStatus.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GENIE.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"

using std::string;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);

//Default options:
int           gOptNevents      = 100;            // n-events to generate
bool          gOptBuildSplines = false;          // spline building option
double        gOptNuEnergy     = 3.0;            // neutrino energy
int           gOptNuPdgCode    = kPdgNuMu;       // neutrino PDG code
int           gOptTgtPdgCode   = 1056026000;     // target PDG code
NtpMCFormat_t gOptNtpFormat    = kNFEventRecord; // ntuple format

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments

  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments

  LOG("test", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("test", pINFO) << "Building splines at init.  = " << gOptBuildSplines;
  LOG("test", pINFO) << "Neutrino energy            = " << gOptNuEnergy;
  LOG("test", pINFO) << "Neutrino PDG code          = " << gOptNuPdgCode;
  LOG("test", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("test", pINFO) << "Output ntuple format       = "
                                    << NtpMCFormat::AsString(gOptNtpFormat);

  //-- create the GENIE high level event generation interface object
  //   for the given initial state

  InitialState init_state(gOptTgtPdgCode, gOptNuPdgCode);

  GENIE genie;
  genie.SetInitialState(init_state);

  //-- load and/or build splines if required
  XSecSplineList * xssl = 0;
  if(gOptBuildSplines) {
     xssl = XSecSplineList::Instance();
     // check whether there is a spline-list XML file to load
     string spllst_load_xmlfile =
             (gSystem->Getenv("GSPLOAD") ? gSystem->Getenv("GSPLOAD") : "");
     LOG("test", pINFO) << "$GSPLOAD env.var = " << spllst_load_xmlfile;

     if(spllst_load_xmlfile.size()>0) {
       LOG("test", pINFO) << "Loading cross section splines from an xml file";
       XmlParserStatus_t status = xssl->LoadFromXml(spllst_load_xmlfile);
       assert(status==kXmlOK);
     }
     // create any spline that is needed but is not loaded
     genie.CreateSplines();
  }

  //-- save the splines if requested
  string spllst_save_xmlfile =
             (gSystem->Getenv("GSPSAVE") ? gSystem->Getenv("GSPSAVE") : "");
  LOG("test", pINFO) << "$GSPSAVE env.var = " << spllst_save_xmlfile;

  if(gOptBuildSplines && spllst_save_xmlfile.size()>0) {
     LOG("test", pINFO) << "Saving cross section splines to an xml file";
     xssl->SaveAsXml(spllst_save_xmlfile);
  }

  //-- in this test just request events for monoenergetic neutrinos
  TLorentzVector nu_p4 (0., 0., gOptNuEnergy, gOptNuEnergy); // px,py,pz,E (GeV)

  //-- create the output ROOT file name;
  ostringstream filename;
  filename << "GNtp" << NtpMCFormat::FilenameTag(gOptNtpFormat) << ".root";

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(gOptNtpFormat);
  ntpw.InitTree(filename.str());

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     // generate a single event
     EventRecord * ev_rec = genie.GenerateEvent(nu_p4);

     // print the event record and the interaction summary
     Interaction & summary = *ev_rec->GetInteraction();

     LOG("test", pINFO) << "Generated Event GHEP Record: " << *ev_rec;
     LOG("test", pINFO) << "Generated Event summary: "     << summary;

     // add event at the output ntuple
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

    if (argv[1][1] == 's') gOptBuildSplines = true;

    int format = 1;
    if (argv[1][1] == 'f') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        format = atoi(argument);
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        format = atoi(argument);
      }
    }
    if(format == 0 || format == 1) gOptNtpFormat = (NtpMCFormat_t)format;

    argc--;
    argv++;
  }
  delete [] argument;
}
//____________________________________________________________________________

