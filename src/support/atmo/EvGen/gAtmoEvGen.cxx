//________________________________________________________________________________________
/*!

\program gevgen_atmo

\brief   A GENIE atmospheric neutrino event generation application.

         *** Syntax :

           gevgen_atmo [-h] 
                       [-r run#] 
                        -f flux
                        -g geometry
                       [-t top_volume_name_at_geom]
                       [-m max_path_lengths_xml_file]
                       [-L length_units_at_geom]
                       [-D density_units_at_geom]
                        -n n_of_events
                        -e min_energy,max_energy
                        [-o output_event_file_prefix]

         *** Options :

           [] Denotes an optional argument

           -h Prints out the syntax and exits

           -r Specifies the MC run number 
              [default: 100000000]

           -f Specifies the input flux files
              The general syntax is: `-f simulation:/path/file.data[neutrino_code],...'
              [Notes] 
               - The `simulation' string can be either `FLUKA' or `BGLRS' (so that
                 input data are binned using the correct FLUKA and BGLRS energy and
                 costheta binning). See comments in 
                 - $GENIE/src/Flux/GFlukaAtmo3DFlux.h
                 - $GENIE/src/Flux/GBartolAtmoFlux.h
                 and follow the links to the FLUKA and BGLRS atmo. flux web pages.
               - The neutrino codes are the PDG ones.
               - The /path/file.data,neutrino_code part of the option can be 
                 repeated multiple times (separated by commas), once for each 
                 flux neutrino species you want to consider, 
                 eg. '-f FLUKA:~/data/sdave_numu07.dat[14],~/data/sdave_nue07.dat[12]'
                 eg. '-f BGLRS:~/data/flux10_271003_z.kam_nue[12]'

           -g Input 'geometry'.
              This option can be used to specify any of:

              1 > A ROOT file containing a ROOT/GEANT geometry description
                  [Note]
                  - This is the standard option for generating events in the
                    nd280, 2km and INGRID detectors.
                  [Examples]
                  - To use the master volume from the ROOT geometry stored
                    in the nd280-geom.root file, type:
                    '-g /some/path/nd280-geom.root'

              2 > A mix of target materials, each with its corresponding weight,
                  typed as a comma-separated list of nuclear pdg codes (in the
                  std PDG2006 convention: 10LZZZAAAI) with the weight fractions
                  in brackets, eg code1[fraction1],code2[fraction2],...
                  If that option is used (no detailed input geometry description)
                  then the interaction vertices are distributed in the detector
                  by the detector MC.
                  [Note]
                  - This is the standard option for generating events in the
                    SuperK detector.
                  [Examples]
                  - To use a target mix of 95% O16 and 5% H type:
                    '-g 1000080160[0.95],1000010010[0.05]'
                  - To use a target which is 100% C12, type:
                    '-g 1000060120'

           -L Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']

           -D Input geometry density units, eg 'g_cm3', 'clhep_def_density_unit',...
              [default: 'g_cm3']

           -t Input 'top volume' for event generation -
              can be used to force event generation in given sub-detector

              [default: the 'master volume' of the input geometry]

              You can also use the -t option to switch generation on/off at
              multiple volumes as, for example, in:
              `-t +Vol1-Vol2+Vol3-Vol4',
              `-t "+Vol1 -Vol2 +Vol3 -Vol4"', 
              `-t -Vol2-Vol4+Vol1+Vol3',
              `-t "-Vol2 -Vol4 +Vol1 +Vol3"'m
              where:
              "+Vol1" and "+Vol3" tells GENIE to `switch on'  Vol1 and Vol3, while
              "-Vol2" and "-Vol4" tells GENIE to `switch off' Vol2 and Vol4.
              If the very first character is a '+', GENIE will neglect all volumes
              except the ones explicitly turned on. Vice versa, if the very first
              character is a `-', GENIE will keep all volumes except the ones
              explicitly turned off (feature contributed by J.Holeczek).

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

           (1) Generate 100k events (run number 999210) in the energy range 1-10 GeV
               for nu_e and nu_mu only, using the sdave_numu07.dat FLUKA flux file for
               nu_mu and the sdave_nue07.dat file for nu_e (files in /data/flx/).
               Use the detector geometry in the /data/geo/SuperK.root file, where the 
               geometry length and density units are m and kgr/m^3. Generate events over 
               the entire geometry volume.

               % gevgen_atmo -r 999210 -n 100000 -e 1,10
                       -f FLUKA:/data/flx/sdave_numu07.dat[14],/data/flx/sdave_nue07.dat[12] 
                       -g /data/geo/SuperK.root -L "m" -D "kg_m3"


           (2) Like above but, instead of generating events in a realistic detector
               geometry, use a simple target mix (88.79% O16 + 11.21% H, i.e. `water')

               % gevgen_atmo -r 999210 -n 100000 -e 1,10
                       -f /data/flux/sdave_numu07.dat[14],/data/flux/sdave_nue07.dat[12] 
                       -g 1000080160[0.8879],1000010010[0.1121]


		... to add more

        
         You can further control the GENIE behaviour by setting its standard 
         environmental variables.
         Please read the GENIE User Manual for more information.

\created August 20, 2010

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Torben Ferber <torben.ferber \at DESY.DE>
         DESY

         Hugh Gallagher <hugh.gallagher \at stfc.ac.uk>
         Tufts University

         Tarak Thakore <tarak \at mailhost.tifr.res.in>
         Tata Institute of Fundamental Research 

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/UnitUtils.h"
#include "Utils/CmdLnArgParser.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GBartolAtmoFlux.h"
#endif

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#include "Geo/GeoUtils.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Geo/PointGeomAnalyzer.h"
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::flux;

void            GetCommandLineArgs (int argc, char ** argv);
void            PrintSyntax        (void);
GFluxI *        GetFlux            (void);
GeomAnalyzerI * GetGeometry        (void);

// User-specified options:
//
Long_t          gOptRunNu;                     // run number
string          gOptFluxSim;                   // flux simulation (FLUKA or BGLRS)
map<int,string> gOptFluxFiles;                 // neutrino pdg code -> flux file map
bool            gOptUsingRootGeom = false;     // using root geom or target mix?
map<int,double> gOptTgtMix;                    // target mix  (tgt pdg -> wght frac) / if not using detailed root geom
string          gOptRootGeom;                  // input ROOT file with realistic detector geometry
string          gOptRootGeomTopVol = "";       // input geometry top event generation volume
double          gOptGeomLUnits = 0;            // input geometry length units
double          gOptGeomDUnits = 0;            // input geometry density units
string          gOptExtMaxPlXml;               // max path lengths XML file for input geometry
int             gOptNev;                       // number of events to generate
double          gOptEvMin;                     // minimum neutrino energy
double          gOptEvMax;                     // maximum neutrino energy
string          gOptEvFilePrefix;              // event file prefix

// Defaults:
//
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // def event tree format
string          kDefOptEvFilePrefix = "gntp";  // def output prefix (override with -o)
string          kDefOptGeomLUnits   = "mm";    // def geom length units (override with -L)
string          kDefOptGeomDUnits   = "g_cm3"; // def geom density units (override with -D)
double          kDefOptEvMin =  0.5;           // min neutrino energy (override with -e)
double          kDefOptEvMax = 50.0;           // max neutrino energy (override with -e)

//________________________________________________________________________________________
int main(int argc, char** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

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
    LOG("gevgen_atmo", pNOTICE) << "Generated event: " << *event;

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
  GeomAnalyzerI * geom_driver = 0;

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__

  if(gOptUsingRootGeom) {
    //
    // *** Using a realistic root-based detector geometry description
    //

    // creating & configuring a root geometry driver
    geometry::ROOTGeomAnalyzer * rgeom =
            new geometry::ROOTGeomAnalyzer(gOptRootGeom);
    rgeom -> SetLengthUnits  (gOptGeomLUnits);
    rgeom -> SetDensityUnits (gOptGeomDUnits);
    rgeom -> SetTopVolName   (gOptRootGeomTopVol);
    // getting the bounding box dimensions along z so as to set the
    // appropriate upstream generation surface for the JPARC flux driver
    TGeoVolume * topvol = rgeom->GetGeometry()->GetTopVolume();
    if(!topvol) {
      LOG("gevgen_atmo", pFATAL) << " ** Null top ROOT geometry volume!";
      gAbortingInErr = true;
      exit(1);
    }
/*
    TGeoShape * bounding_box = topvol->GetShape();
    bounding_box->GetAxisRange(3, zmin, zmax);
    zmin *= rgeom->LengthUnits();
    zmax *= rgeom->LengthUnits();
*/
    // switch on/off volumes as requested
    if ( (gOptRootGeomTopVol[0] == '+') || (gOptRootGeomTopVol[0] == '-') ) {
      bool exhaust = (*gOptRootGeomTopVol.c_str() == '+');
      utils::geometry::RecursiveExhaust(topvol, gOptRootGeomTopVol, exhaust);
    }

    // casting to the GENIE geometry driver interface
    geom_driver = dynamic_cast<GeomAnalyzerI *> (rgeom);
  }
  else {
    //       
    // *** Using a 'point' geometry with the specified target mix
    // *** ( = a list of targets with their corresponding weight fraction)
    //
  
    // creating & configuring a point geometry driver
    geometry::PointGeomAnalyzer * pgeom =
              new geometry::PointGeomAnalyzer(gOptTgtMix);
    // casting to the GENIE geometry driver interface
    geom_driver = dynamic_cast<GeomAnalyzerI *> (pgeom);
  }

#else
  LOG("gevgen_atmo", pFATAL) << "You need to enable the GENIE geometry drivers first!";
  LOG("gevgen_atmo", pFATAL) << "Use --enable-geom-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return geom_driver;
}
//________________________________________________________________________________________
GFluxI* GetFlux(void)
{
  GFluxI * flux_driver = 0;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

  // Instantiate appropriate concrete flux driver
  GAtmoFlux * atmo_flux_driver = 0;
  if(gOptFluxSim == "FLUKA") {
     GFlukaAtmo3DFlux * fluka_flux = new GFlukaAtmo3DFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(fluka_flux);
  } else
  if(gOptFluxSim == "BGLRS") {
     GBartolAtmoFlux * bartol_flux = new GBartolAtmoFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(bartol_flux);
  } else {
     LOG("gevgen_atmo", pFATAL) << "Uknonwn flux simulation: " << gOptFluxSim;
     gAbortingInErr = true;
     exit(1);
  }
  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set min/max energy:
  atmo_flux_driver->ForceMinEnergy (gOptEvMin * units::GeV);
  atmo_flux_driver->ForceMaxEnergy (gOptEvMax * units::GeV);    
  // set flux files:
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
    int neutrino_code = file_iter->first;
    string filename   = file_iter->second;
    atmo_flux_driver->SetFluxFile(neutrino_code, filename);
  }
  atmo_flux_driver->LoadFluxData();
  atmo_flux_driver->SetRadii(1, 1);
  // Cast to GFluxI, the generic flux driver interface 
  flux_driver = dynamic_cast<GFluxI *>(atmo_flux_driver);

#else
  LOG("gevgen_atmo", pFATAL) << "You need to enable the GENIE flux drivers first!";
  LOG("gevgen_atmo", pFATAL) << "Use --enable-flux-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return flux_driver;
}
//________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Get the command line arguments

  LOG("gevgen_atmo", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  //
  // *** run number
  //
  if( parser.OptionExists('r') ) {
    LOG("gevgen_atmo", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_atmo", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 100000000;
  } //-r

  //
  // *** exposure
  // 

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_atmo", pDEBUG) 
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_atmo", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  } //-n

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_atmo", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_atmo", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  //
  // *** neutrino energy range
  //
  if( parser.OptionExists('e') ) {
    LOG("gevgen_atmo", pINFO) << "Reading neutrino energy range";
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
      LOG("gevgen_atmo", pFATAL)
        << "Invalid energy range. Use `-e emin,emax', eg `-e 0.5,100.";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  } else {
     LOG("gevgen_atmo", pNOTICE)
        << "No -e option. Using default energy range";
     gOptEvMin = kDefOptEvMax;
     gOptEvMax = kDefOptEvMax;
  }

  //
  // *** flux files
  //
  // syntax: 
  // simulation:/path/file.data[neutrino_code],/path/file.data[neutrino_code],...
  //
  if( parser.OptionExists('f') ) {
    LOG("gevgen_atmo", pDEBUG) << "Getting input flux files";
    string flux = parser.ArgAsString('f');

    // get flux simulation info (FLUKA,BGLRS) so as to instantiate the
    // appropriate flux driver
    string::size_type jsimend = flux.find_first_of(":",0);
    if(jsimend==string::npos) {
       LOG("gevgen_atmo", pFATAL) 
           << "You need to specify the flux file source"; 
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }
    gOptFluxSim = flux.substr(0,jsimend);
    for(string::size_type i=0; i<gOptFluxSim.size(); i++) {
       gOptFluxSim[i] = toupper(gOptFluxSim[i]);
    }
    if((gOptFluxSim != "FLUKA") && (gOptFluxSim != "BGLRS")) {
        LOG("gevgen_atmo", pFATAL) 
             << "The flux file source needs to be one of <FLUKA,BGLRS>"; 
        PrintSyntax();
        gAbortingInErr = true;
        exit(1);
    }
    // now get the list of input files and the corresponding neutrino codes.
    flux.erase(0,jsimend+1); 
    vector<string> fluxv = utils::str::Split(flux,",");      
    vector<string>::const_iterator fluxiter = fluxv.begin();
    for( ; fluxiter != fluxv.end(); ++fluxiter) {
       string filename_and_pdg = *fluxiter;
       string::size_type open_bracket  = filename_and_pdg.find("[");
       string::size_type close_bracket = filename_and_pdg.find("]");
       if (open_bracket ==string::npos || 
           close_bracket==string::npos) 
       {
           LOG("gevgen_atmo", pFATAL) 
              << "You made an error in specifying the flux info"; 
           PrintSyntax();
           gAbortingInErr = true;
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
    if(gOptFluxFiles.size() == 0) {
       LOG("gevgen_atmo", pFATAL) 
          << "You must specify at least one flux file!"; 
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }

  } else {
    LOG("gevgen_atmo", pFATAL) << "No flux info was specified! Use the -f option.";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  //
  // *** geometry
  //        
  string geom = "";
  string lunits, dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_atmo", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');
     
    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
        utils::system::FileExists(geom.c_str());
    if (accessible_geom_file) {
      gOptRootGeom      = geom; 
      gOptUsingRootGeom = true;
    }
  } else {  
      LOG("gevgen_atmo", pFATAL)
        << "No geometry option specified - Exiting";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units
    
     // legth units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_atmo", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_atmo", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // density units:
     if( parser.OptionExists('D') ) {
        LOG("gevgen_atmo", pDEBUG)
           << "Checking for input geometry density units";
        dunits = parser.ArgAsString('D');
     } else {
        LOG("gevgen_atmo", pDEBUG) << "Using default geometry density units";
        dunits = kDefOptGeomDUnits;
     } // -D
     gOptGeomLUnits = genie::utils::units::UnitFromString(lunits);
     gOptGeomDUnits = genie::utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_atmo", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_atmo", pDEBUG) << "Using the <master volume>";
     } // -t

     // check whether an XML file with the maximum (density weighted)
     // path lengths for each detector material is specified -
     // otherwise will compute the max path lengths at job init
     if( parser.OptionExists('m') ) {
        LOG("gevgen_atmo", pDEBUG)
              << "Checking for maximum path lengths XML file";
        gOptExtMaxPlXml = parser.ArgAsString('m');
     } else {
        LOG("gevgen_atmo", pDEBUG)
           << "Will compute the maximum path lengths at job init";
        gOptExtMaxPlXml = "";
     } // -m
  } // using root geom?

  else {
    // User has specified a target mix.
    // Decode the list of target pdf codes & their corresponding weight fraction
    // (specified as 'pdg_code_1[fraction_1],pdg_code_2[fraction_2],...')
    // See documentation on top section of this file.
    //
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(geom,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } else {
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
         string tgt_with_wgt = *tgtmix_iter;
         string::size_type open_bracket  = tgt_with_wgt.find("[");
         string::size_type close_bracket = tgt_with_wgt.find("]");
         if (open_bracket ==string::npos ||
             close_bracket==string::npos)
         {
             LOG("gevgen_atmo", pFATAL)
                << "You made an error in specifying the target mix";
             PrintSyntax();
             gAbortingInErr = true;
             exit(1);
         }
         string::size_type ibeg = 0;
         string::size_type iend = open_bracket;
         string::size_type jbeg = open_bracket+1;
         string::size_type jend = close_bracket;
         int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend-ibeg).c_str());
         double wgt = atof(tgt_with_wgt.substr(jbeg,jend-jbeg).c_str());
         LOG("gevgen_atmo", pDEBUG)
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));

      }// tgtmix_iter
    } // >1 materials in mix
  } // using tgt mix?
    
  //
  // print-out summary
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();

  ostringstream gminfo;
  if (gOptUsingRootGeom) {
    gminfo << "Using ROOT geometry - file: " << gOptRootGeom
           << ", top volume: "
           << ((gOptRootGeomTopVol.size()==0) ? "<master volume>" : gOptRootGeomTopVol)
           << ", max{PL} file: "
           << ((gOptExtMaxPlXml.size()==0) ? "<none>" : gOptExtMaxPlXml)
           << ", length  units: " << lunits
           << ", density units: " << dunits;
  } else {
    gminfo << "Using target mix - ";
    map<int,double>::const_iterator iter;
    for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
          int    pdg_code = iter->first;
          double wgt      = iter->second;
          TParticlePDG * p = pdglib->Find(pdg_code);
          if(p) {
            string name = p->GetName();
            gminfo << "(" << name << ") -> " << 100*wgt << "% / ";
          }//p?
    }
  }
     
  ostringstream fluxinfo;
  fluxinfo << "Using " << gOptFluxSim << " flux files: ";
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
     int neutrino_code = file_iter->first;
     string filename   = file_iter->second;
     TParticlePDG * p = pdglib->Find(neutrino_code);
     if(p) {
        string name = p->GetName();
        fluxinfo << "(" << name << ") -> " << filename << " / ";
     }
  }

  LOG("gevgen_atmo", pNOTICE) 
     << "\n"
     << "\n ****** MC Job (" << gOptRunNu << ") Settings ****** "
     << "\n @@ Geometry"
     << "\n\t" << gminfo.str()
     << "\n @@ Flux"
     << "\n\t" << fluxinfo.str()
     << "\n @@ Exposure" 
     << "\n\t Number of events = " << gOptNev
     << "\n @@ Cuts"
     << "\n\t Using energy range = (" << gOptEvMin << " GeV, " << gOptEvMax << " GeV)"
     << "\n\n";
}
//________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_atmo", pFATAL) 
   << "\n **Syntax**"
   << "\n gevgen_atmo [-h]"
   << "\n           [-r run#]"
   << "\n            -f simulation:flux_file[neutrino_code],..."
   << "\n            -g geometry"
   << "\n           [-t top_volume_name_at_geom]"
   << "\n           [-m max_path_lengths_xml_file]"
   << "\n           [-L length_units_at_geom]"
   << "\n           [-D density_units_at_geom]"
   << "\n            -n n_of_events"
   << "\n            -e min_energy,max_energy"
   << "\n           [-o output_event_file_prefix]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n";
}
//________________________________________________________________________________________

