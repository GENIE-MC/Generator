//____________________________________________________________________________
/*!

\program gmkspl

\brief   GENIE utility program building XML cross section splines that can be
         loaded into GENIE to speed-up event generation

         Syntax :
           gmkspl -p nupdg -t tgtpdg [-f output_xml_file]

         Options :
           -p a comma separated list of nu PDG code
           -t a comma separated list of tgt PDG codes (format: 1aaazzz000)
           -f output XML filename

         Example:
           gmkspl -p 14,-14 -t 1056026000

           will build cross section splines for muon neutrinos (pdg = 14)
           and muon anti-neutrinos (pgc = -14) on Iron (A=56,Z=26).

         Other control options:

         You can further control the program behaviour by setting the GEVGL
         and GMSGCONF environmental variables.

           - Set the GEVGL environmental variable to contol the list of event
             generator objects that will get loaded into the event generation
             driver (see $GENIE/src/stdapp/gEvGen.cxx). This program will only
             build splines for the processes that can be simulated by the event
             generators you plan to load.

           - You can set the GMSGCONF env.variable to point to a messenger
             XML configuration file (following the syntax of the default one
             that can be found in $GENIE/config/messenger.xml) and modify the
             verbosity of GENIE output. Both the default and your messenger
             configuration will be read but yours will take precedence in
             case of clashing priority for the same stream.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 27, 2005
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <vector>

#include <TSystem.h>

#include "EVGDrivers/GEVGDriver.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/XSecSplineList.h"

using std::string;
using std::vector;

using namespace genie;

void GetCommandLineArgs(int argc, char ** argv);

//Default options:
string gOptNuPdgCodeList  = "";  // comma separated list of neutrino PDG codes
string gOptTgtPdgCodeList = "";  // comma separated list of target PDG codes
string gOptXMLFilename = "xsec_splines.xml";

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments
  LOG("gmkspl", pINFO) << "Neutrino PDG codes = " << gOptNuPdgCodeList;
  LOG("gmkspl", pINFO) << "Target PDG codes   = " << gOptTgtPdgCodeList;
  LOG("gmkspl", pINFO) << "Output XML file    = " << gOptXMLFilename;

  //-- split the coma-separated lists
  vector<string> nuvec  = string_utils::Split(gOptNuPdgCodeList,  ",");
  vector<string> tgtvec = string_utils::Split(gOptTgtPdgCodeList, ",");

  if(nuvec.size() == 0 || tgtvec.size() == 0) {

    if(nuvec.size() == 0) {
         LOG("gmkspl", pFATAL) << "Empty neutrino PDG code list";
    }
    if(tgtvec.size() == 0) {
         LOG("gmkspl", pFATAL) << "Empty target PDG code list";
    }
    LOG("gmkspl", pINFO)
      << "\n\n Syntax: gmkspl -p nupdg -t tgtpdg [-f output_xml_file] \n\n";
  }

  //-- loop over all possible input init states and ask the GEVGDriver
  //   to build splines for all the interactions that its loaded list
  //   of event generators can generate.

  vector<string>::const_iterator nuiter;
  vector<string>::const_iterator tgtiter;

  for(nuiter = nuvec.begin(); nuiter != nuvec.end(); ++nuiter) {
    for(tgtiter = tgtvec.begin(); tgtiter != tgtvec.end(); ++tgtiter) {

      int nupdgc  = atoi(nuiter->c_str());
      int tgtpdgc = atoi(tgtiter->c_str());

      InitialState init_state(tgtpdgc, nupdgc);

      GEVGDriver driver;
      driver.SetInitialState(init_state);

      driver.CreateSplines();
    }
  }

  //-- get the populated cross section spline list and save it at the
  //   requested XML file

  XSecSplineList * spline_list = XSecSplineList::Instance();

  spline_list->SaveAsXml(gOptXMLFilename);

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  char * argument = new char[1024];

  while( argc>1 && (argv[1][0] == '-'))
  {
    if (argv[1][1] == 'p') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptNuPdgCodeList = argument;
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptNuPdgCodeList = argument;
      }
    }
    if (argv[1][1] == 't') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptTgtPdgCodeList = argument;
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptTgtPdgCodeList = argument;
      }
    }
    if (argv[1][1] == 'f') {
      if (strlen(&argv[1][2]) ) {
        strcpy(argument,&argv[1][2]);
        gOptXMLFilename = argument;
      } else if( (argc>2) && (argv[2][0] != '-') ) {
        argc--;
        argv++;
        strcpy(argument,&argv[1][0]);
        gOptXMLFilename = argument;
      }
    }
    argc--;
    argv++;
  }
  delete [] argument;
}
//____________________________________________________________________________

