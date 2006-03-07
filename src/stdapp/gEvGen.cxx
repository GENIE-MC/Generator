//____________________________________________________________________________
/*!

\program gevgen

\brief   Example program driving GENIE event generation modules

         Syntax :
           gevgen [-n nev] [-s] -e E -p nupdg -t tgtpdg [-f format] [-r run#]

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate
           -s turns on cross section spline building at job initialization
           -e specifies the neutrino energy
	      (if what follows the -e option is a comma separated pair of
               values it will be interpreted as an energy range and neutrino
               energies will be generated uniformly in this range).
              Note that this is the energy of an "interacting" neutrino, not 
              the energy of a "flux" neutrino -- the  neutrino is forced to 
              interact anyway (see GMCJDriver for how to input fluxes)
           -p specifies the neutrino PDG code
           -t specifies the target PDG code (std format: 1aaazzz000)
           -f specifies the output TTree format. If set to 0 will create a
              single-branch TTree with NtpMCPlainRecord objects in its leaves,
              while if set to 1 it will have NtpMCEventRecord objects in
              its leaves (see the Ntuple package for descriptions of the ntuple
              records and their intended usage). Default options is 1.
           -r specifies the MC run number

         Example:
           gevgen -n 300 -s -e 6.5 -p 14 -t 1056026000 -r 10000

           will build cross section splines during the initialization step,
           and will generate 300 events of muon neutrinos (pdg = 14) on Iron
           (A=56,Z=26) at E = 6.5 GeV.

         Other control options:

         You can further control the program behaviour by setting the GEVGL,
         GSPLOAD, GSPSAVE, GHEPPRINTLEVEL, GSEED environmental variables.

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

           - You can set the GHEPPRINTLEVEL to control the GHEP print options

           - You can set the GSEED env variable to define the random number
             seed number (default seed number in (src/Conventions/Controls.h)

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
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "Conventions/XmlParserStatus.h"
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

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 100;            // n-events to generate
NtpMCFormat_t kDefOptNtpFormat = kNFEventRecord; // ntuple format
Long_t        kDefOptRunNu     = 0;              // default run number

//User-specified options:
int           gOptNevents;      // n-events to generate
bool          gOptBuildSplines; // spline building option
double        gOptMinNuEnergy;  // min neutrino energy 
double        gOptNuEnergyRange;// max-min neutrino energy 
int           gOptNuPdgCode;    // neutrino PDG code
int           gOptTgtPdgCode;   // target PDG code
NtpMCFormat_t gOptNtpFormat;    // ntuple format
Long_t        gOptRunNu;        // run number

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);
  
  //-- print the options you got from command line arguments

  string fmts = NtpMCFormat::AsString(gOptNtpFormat);

  LOG("gevgen", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("gevgen", pINFO) << "Building splines at init.  = " << gOptBuildSplines; 
  LOG("gevgen", pINFO) << "Neutrino PDG code          = " << gOptNuPdgCode;
  LOG("gevgen", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("gevgen", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("gevgen", pINFO) << "Output ntuple format       = " << fmts;
  if(gOptNuEnergyRange>0) {
    LOG("gevgen", pINFO) << "Neutrino energy            = [" 
        << gOptMinNuEnergy << ", " << gOptMinNuEnergy+gOptNuEnergyRange << "]";
  } else {
    LOG("gevgen", pINFO) << "Neutrino energy            = " << gOptMinNuEnergy;
  }

  //-- create the GENIE high level event generation interface object
  //   for the given initial state

  InitialState init_state(gOptTgtPdgCode, gOptNuPdgCode);

  GEVGDriver driver;
  driver.Configure(init_state);

  //-- Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  //   if the env. var. has been set)
  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->AutoLoad();

  //-- If splines are used, then create any spline that is needed but is not
  //   already loaded (can built them all here if no spline at all is loaded)
  if(gOptBuildSplines) driver.CreateSplines();

  //-- this driver produces events for monoenergetic neutrinos
  TLorentzVector nu_p4(0.,0.,gOptMinNuEnergy,gOptMinNuEnergy); // px,py,pz,E (GeV)

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(gOptNtpFormat, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get a random number generator
  RandomGen * r = RandomGen::Instance();

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     // generate neutrino energy (if an energy range was defined)
     if(gOptNuEnergyRange>0) {
       double Ev = gOptMinNuEnergy + gOptNuEnergyRange * r->Random1().Rndm();
       nu_p4.SetPxPyPzE(0,0,Ev,Ev);
     }

     // generate a single event
     EventRecord * event = driver.GenerateEvent(nu_p4);

     LOG("gevgen", pINFO) << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple
     ntpw.AddEventRecord(ievent, event);

     // refresh the mc job monitor
     mcjmonitor.Update(ievent,event);

     ievent++;
     delete event;
  }

  //-- save the generated MC events
  ntpw.Save();

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen", pNOTICE) << "Parsing commad line arguments";

  //-- Optional arguments

  //number of events:
  try {
    LOG("gevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  //run number:
  try {
    LOG("gevgen", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  //output ntuple format
  int format = 1;
  try {
    LOG("gevgen", pINFO) << "Reading requested output ntuple format";
    format = genie::utils::clap::CmdLineArgAsInt(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO) << "Unspecified tree format - Using default";
    }
  }
  if(format == 0 || format == 1) gOptNtpFormat = (NtpMCFormat_t)format;

  //spline building option
  gOptBuildSplines = genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  //-- Required arguments

  //neutrino energy:
  try {
    LOG("gevgen", pINFO) << "Reading neutrino energy";
    string nue = genie::utils::clap::CmdLineArgAsString(argc,argv,'e');

    // is it just a value or a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);   
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>0);
       gOptMinNuEnergy   = emin;
       gOptNuEnergyRange = emax-emin;
    } else {
       gOptMinNuEnergy   = atof(nue.c_str());
       gOptNuEnergyRange = -1;
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified neutrino energy - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //neutrino PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //target PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevgen [-n nev] [-s] -e energy -p nupdg -t tgtpdg [-f format]\n";
}
//____________________________________________________________________________
