//________________________________________________________________________________________
/*!

\program gLBNEAtmosEvGen

\brief   A GENIE event generation driver 'customized' for LBNE, modified from INO

         *** Syntax :

           gLBNEAtmosEvGen [-h] 
                           [-r run#] 
                           -f flux
                           -n n_of_events
                           -e min_energy,max_energy
                           [-o output_event_file_prefix]

         *** Options :

           [] Denotes an optional argument

           -h Prints out the syntax and exits

           -r Specifies the MC run number 
              [default: 0]

           -f Specifies the input flux files
              The general syntax is: `-f /path/file.data[neutrino_code],...'
              [Notes] 
               - The neutrino codes are the PDG ones.
               - The /path/file.data,neutrino_code part of the option can be 
                 repeated multiple times (separated by commas), once for each 
                 flux neutrino species you want to consider, 
                 eg. '-f ~/data/sdave_numu07.dat[14],~/data/sdave_nue07.dat[12]'

           -n Specifies how many events to generate.

           -e Specifies the neutrino energy in GeV. 
              Must be a comma-separated pair of numbers, eg `-e 0.3,70'
              [default: 0.5,50]

           -o Sets the prefix of the output event file. 
              The output filename is built as: 
              [prefix].[run_number].[event_tree_format].[file_format]
              The default output filename is: 
              gntp.[run_number].ghep.root
              This cmd line arguments lets you override 'gntp'

         *** Examples:


           (1) Generate 10k events (run number 999210) in the energy range 1-10 GeV
               for nu_e and nu_mu only, using the sdave_numu07.dat FLUKA flux file
               for nu_mu's and the sdave_nue07.dat FLUKA flux file for nu_e's:

               % gLBNEAtmosEvGen -r 999210 -n 10000 -e 1,10
                       -f /data/flux/sdave_numu07.dat[14],/data/flux/sdave_nue07.dat[12] 


		... to add more
        
         You can further control the GENIE behaviour by setting its standard 
         environmental variables.
         Please read the GENIE User Manual for more information.

\created July 16, 2010

\author  Tarak Thakore <tarak \at mailhost.tifr.res.in>
         TATA INSTITUTE OF FUNDAMENTAL RESEARCH 

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Hugh Gallagher <hugh.gallagher \at stfc.ac.uk>
         Tufts University

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GBartolAtmoFlux.h"
#include "Geo/PointGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "PDG/PDGCodes.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using namespace std;
using namespace genie;
using namespace genie::flux;

void            GetCommandLineArgs (int argc, char ** argv);
void            PrintSyntax        (void);
GFluxI *        GetFlux            (void);
GeomAnalyzerI * GetGeometry        (void);

// User-specified options:
//
Long_t          gOptRunNu;                     // run number
int             gOptNev;                       // number of events to generate
double          gOptEvMin;                     // minimum neutrino energy
double          gOptEvMax;                     // maximum neutrino energy
string          gOptEvFilePrefix;              // event file prefix
map<int,string> gOptFluxFiles;                 // neutrino pdg code -> flux file map

// Defaults:
//
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // deflt event tree format
string          kDefOptEvFilePrefix = "gntp";  // deflt output prefix (override with -o)
double          kDefOptEvMin =  0.5;           // min neutrino energy (override with -e)
double          kDefOptEvMax = 50.0;           // max neutrino energy (override with -e)

//________________________________________________________________________________________
int main(int argc, char** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Set seed
  RandomGen::Instance()->SetSeed(gOptRunNu);

  // Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  // if the env. var. has been set)
  XSecSplineList::Instance()->AutoLoad();

  // get flux driver
  GFluxI * flux_driver = GetFlux();

  // get geometry driver
  GeomAnalyzerI * geom_driver = GetGeometry();

  // create the GENIE monte carlo job driver
  GMCJDriver* mcj_driver = new GMCJDriver;
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  mcj_driver->ForceSingleProbScale();

  // initialize an ntuple writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);

  // event loop
  for(int iev = 0; iev < gOptNev; iev++) {

    // generate next event
    EventRecord* event = mcj_driver->GenerateEvent();

    // set weight (if using a weighted flux)
    //event->SetWeight(event->Weight()*flux_driver->Weight());

    // print-out
    LOG("gLBNEAtmosEvGen", pNOTICE) 
        << "Generated event: " << *event;

    // save the event, refresh the mc job monitor 
    ntpw.AddEventRecord(iev, event);
    mcjmonitor.Update(iev,event);

    // clean-up
    delete event;
  }

  // save the event file
  ntpw.Save();

  // clean-up
  delete geom_driver;
  delete flux_driver;
  delete mcj_driver;

  return 0;
}
//________________________________________________________________________________________
GeomAnalyzerI* GetGeometry(void)
{
  // create a trivial point geometry
  int targetCode = kPdgTgtC12;
  GeomAnalyzerI* geom_driver = new geometry::PointGeomAnalyzer(targetCode);
  return geom_driver;
}
//________________________________________________________________________________________
GFluxI* GetFlux(void)
{
  GFlukaAtmo3DFlux * flux = new GFlukaAtmo3DFlux;

  // set min/max energy
  flux->ForceMinEnergy (gOptEvMin /* GeV  */);
  flux->ForceMaxEnergy (gOptEvMax /* GeV  */);    

  // set flux files
  assert(gOptFluxFiles.size() != 0);

  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
     int neutrino_code = file_iter->first;
     string filename   = file_iter->second;
     flux->SetFluxFile(neutrino_code, filename);
  }
  flux->LoadFluxData();
  flux->SetRadii(1, 1);

  GFluxI * flux_base = dynamic_cast<GFluxI *>(flux);
  return flux_base;
}
//________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  //
  // >>> get the command line arguments
  //

  LOG("gLBNEAtmosEvGen", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // run number:
  if( parser.OptionExists('r') ) {
    LOG("gLBNEAtmosEvGen", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gLBNEAtmosEvGen", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 0;
  } //-r

  // number of events:
  if( parser.OptionExists('n') ) {
    LOG("gLBNEAtmosEvGen", pDEBUG) 
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gLBNEAtmosEvGen", pERROR)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // event file prefix:
  if( parser.OptionExists('o') ) {
    LOG("gLBNEAtmosEvGen", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gLBNEAtmosEvGen", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // neutrino energy range:
  if( parser.OptionExists('e') ) {
    LOG("gLBNEAtmosEvGen", pINFO) << "Reading neutrino energy range";
    string nue = parser.ArgAsString('e');

    // must be a comma separated set of values
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0.);
       gOptEvMin = emin;
       gOptEvMax = emax;
    } else {
      LOG("gLBNEAtmosEvGen", pERROR)
        << "Invalid energy range. Use `-e emin,emax', eg `-e 0.5,100.";
      PrintSyntax();
      exit(0);
    }
  } else {
     LOG("gLBNEAtmosEvGen", pNOTICE)
        << "No -e option. Using default energy range";
     gOptEvMin = kDefOptEvMax;
     gOptEvMax = kDefOptEvMax;
  }

  // flux files:
  if( parser.OptionExists('f') ) {
    LOG("gLBNEAtmosEvGen", pDEBUG) << "Getting input flux files";
    string flux = parser.ArgAsString('f');

    vector<string> fluxv = utils::str::Split(flux,",");      
    vector<string>::const_iterator fluxiter = fluxv.begin();
    for( ; fluxiter != fluxv.end(); ++fluxiter) {
       string filename_and_pdg = *fluxiter;
       string::size_type open_bracket  = filename_and_pdg.find("[");
       string::size_type close_bracket = filename_and_pdg.find("]");
       if (open_bracket ==string::npos || 
           close_bracket==string::npos) 
       {
           LOG("gLBNEAtmosEvGen", pERROR) 
              << "You made an error in specifying the flux info"; 
           PrintSyntax();
           exit(1);
       }
       string::size_type ibeg = 0;
       string::size_type iend = open_bracket;
       string::size_type jbeg = open_bracket+1;
       string::size_type jend = close_bracket;
       string flux_filename   = filename_and_pdg.substr(ibeg,iend-ibeg);
       string neutrino_pdg    = filename_and_pdg.substr(jbeg,jend-jbeg);
       gOptFluxFiles.insert( 
          map<int,string>::value_type(atoi(neutrino_pdg.c_str()), flux_filename));
    }

  } else {
    LOG("gLBNEAtmosEvGen", pFATAL) << "No flux info was specified! Use the -f option.";
    PrintSyntax();
    exit(1);
  }

  // print out
  LOG("gLBNEAtmosEvGen", pNOTICE) << "Job configuration:";
  LOG("gLBNEAtmosEvGen", pNOTICE) << "Number of events = " << gOptNev;
  LOG("gLBNEAtmosEvGen", pNOTICE) << "Emin = " << gOptEvMin;
  LOG("gLBNEAtmosEvGen", pNOTICE) << "Emax = " << gOptEvMax;
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
     int neutrino_code = file_iter->first;
     string filename   = file_iter->second;
     LOG("gLBNEAtmosEvGen", pNOTICE) 
       << "Neutrino = " << neutrino_code << " --> Flux file = " << filename;
  }

  LOG("gLBNEAtmosEvGen", pNOTICE) << "Done reading command line arguments";
}
//________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gLBNEAtmosEvGen", pFATAL) 
   << "\n **Syntax**"
   << "\n gLBNEevgen_atmo [-h] "
   << "\n           [-r run#]"
   << "\n            -f flux_file[neutrino_code],..."
   << "\n            -n n_of_events"
   << "\n            -e min_energy,max_energy"
   << "\n           [-o output_event_file_prefix]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n";
}
//________________________________________________________________________________________

