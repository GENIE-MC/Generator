//____________________________________________________________________________
/*!

\program testMCJobDriver

\brief   Simple program to drive the GMCJobDriver and generate events for an
         input neutrino flux and an input detector geometry

         Syntax :
           testMCJobDriver -f filename [-u units] [-r run] [-n nevents] [-s]
                          [-m filename] [-d dx,dy,dz] [-b x,y,z] [-t size]

         Options :
           -f  a ROOT file containing a ROOT/GEANT geometry description
           -u  geometry length units [default: meter]
           -n  number of events [default: 10]
           -r  run number [default: 0]
           -m  XML file with max path lengths for geometry (if not set then
               the max path lengths for the input geometry will be computed
               at the mc job driver configuration and it might take a while
               for complex ROOT/GEANT geometries)
           -s  if set turns on using xsec splines (if set splines will be
               loaded from the $GSPLOAD XML file and all needed but missing
               splines will be computed)
           -d  beam direction (dx,dy,dz) [default: 1,0,0 - along x]
           -b  beam "spot" (x,y,z) - defines a surface perpendicular to beam
               direction where flux neutrinos will be generated  (in m)
               [default: -10m, 0m, 0m]
           -t  radius of the area, on the surface perpendicular to the beam
               direction, in which flux neutrinos would be generated (in m)
               [default: 5 m]

         Notes:
           - An argument in []'s is optional
           - The defaults for the flux driver match the example detector setup
             in the $GENIE/src/test/data directory

         Example :
           testMCJobDriver -f ~/mysim/inputs/geometry.root -u cm -r 101 -n 1000
                           -s -m ~/mysim/inputs/maxpl.xml -d 0,0,1 -b 0,0,-100
                           -t 0.20

         Also see $GENIE/src/stdapp/gEvGen.cxx for a list of environmental
         variables that can affect the program behaviour (eg which event
         generators to load, how to specify the xsec splines XML file,...)
         For more options for the used flux driver and more flux drivers see
         $GENIE/src/FluxDrivers/

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 22, 2005
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "PDG/PDGCodes.h"
#include "Utils/XSecSplineList.h"
#include "Utils/UnitUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using std::vector;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

//prototypes:
void            PrintSyntax        (void);
void            GetCommandLineArgs (int argc, char ** argv);
GFluxI *        GetFluxDriver      (void);
GeomAnalyzerI * GetGeometryDriver  (void);

//command line options:
Long_t   gOptRunNu;         // run number
bool     gOptBuildSplines;  // spline building option
int      gOptNevents;       // number of events to generate
string   gOptRootGeom;      // detector geometry ROOT file
string   gOptGeomUnits;     // detector geometry units
string   gOptExtMaxPlXml;   // external max path lengths XML file
TVector3 gOptBeamDirection; // beam direction
TVector3 gOptBeamSpot;      // beam spot
double   gOptBeamRt;        // beam transverse size (radius)

// defaults for optional command line arguments
Long_t   kDefOptRunNu     = 0;
string   kDefOptGeomUnits = "m";
int      kDefOptNevents   = 10;
double   kDefOptBeamRt    = 5;
TVector3 kDefOptBeamDirection(1,0,0);
TVector3 kDefOptBeamSpot(-10,0,0);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- Parse command line arguments
  GetCommandLineArgs(argc, argv);

  //-- Create/configure a flux driver
  GFluxI * flux = GetFluxDriver();

  //-- Create/configure a geometry driver
  GeomAnalyzerI * geom = GetGeometryDriver();

  //-- Create/configure the GENIE MC-job driver
  LOG("Main", pINFO) << "Creating the GENIE MC Job driver";
  GMCJDriver mcj;

  LOG("Main", pINFO) << "Loading flux & geometry drivers";
  mcj.UseFluxDriver(flux);
  mcj.UseGeomAnalyzer(geom);
  mcj.UseMaxPathLengths(gOptExtMaxPlXml);

  LOG("Main", pINFO) << "Configuring the GENIE MC Job driver";
  mcj.Configure();

  //-- If this job uses cross section splines, build all splines that
  //   are needed and have not already loaded from an XML file via
  //   XSecSplineList::AutoLoad()
  if(gOptBuildSplines) mcj.UseSplines();

  //-- Initialize an Ntuple Writer
  NtpWriter ntpw(kNFEventRecord);
  ntpw.Initialize("mcjobdriver.root");

  //-- Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- Event loop
  int i=0;
  while (i<gOptNevents) {

     // generate next event
     EventRecord * event = mcj.GenerateEvent();

     // update job monitor
     mcjmonitor.Update(i,event);

     // print event & add it at the output event tree
     if(event) {
       LOG("Main", pINFO) << "\n\n***EVENT NU: " << i;
       LOG("Main", pINFO) << *event;

       ntpw.AddEventRecord(i++, event);
       delete event;

     } else {
       LOG("Main", pINFO)
         << "\n** Got a null interaction."
         << "If recursive mode isn't allowed then an error occurred";
     }
  }

  //-- Save the generated MC events
  ntpw.Save();

  //-- Clean up and exit
  delete flux;
  delete geom;

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
GFluxI * GetFluxDriver(void)
{
// Create/configure a flux driver

  LOG("Main", pINFO)
            << "Creating/configuring the GCylindTH1Flux flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  TF1 * f1 = new TF1("f1","1./x+2.",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu spectrum", 20,0.5,5);
  spectrum1->FillRandom("f1",10000);

  flux -> SetNuDirection      (gOptBeamDirection);
  flux -> SetBeamSpot         (gOptBeamSpot);
  flux -> SetTransverseRadius (gOptBeamRt);
  flux -> AddEnergySpectrum   (kPdgNuMu, spectrum1);

  GFluxI * fluxb = dynamic_cast<GFluxI *> (flux);

  delete f1;

  return fluxb;
}
//___________________________________________________________________
GeomAnalyzerI * GetGeometryDriver(void)
{
// Create/configure a geometry driver

  LOG("Main", pINFO) << "Creating/configuring the ROOT geom. driver";

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptRootGeom);
  geom->SetUnits(genie::utils::units::UnitFromString(gOptGeomUnits));

  GeomAnalyzerI * geomb = dynamic_cast<GeomAnalyzerI *> (geom);

  return geomb;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Parse the command line arguments

  //geometry file:
  try {
    LOG("Main", pINFO) << "Getting input geometry file";
    gOptRootGeom =
              genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pFATAL) << "No geometry file was specified - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
  //geometry units:
  try {
    LOG("Main", pINFO) << "Getting input geometry units";
    gOptGeomUnits =
              genie::utils::clap::CmdLineArgAsString(argc,argv,'u');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pNOTICE) << "Using default geometry units";
      gOptGeomUnits = kDefOptGeomUnits;
    }
  }
  //external, maximum path lengths:
  try {
    LOG("Main", pINFO) << "Getting XML file with max path length list";
    gOptExtMaxPlXml =
              genie::utils::clap::CmdLineArgAsString(argc,argv,'m');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pFATAL)
        << "No input max path lengths. They will be computed at jon init";
     gOptExtMaxPlXml = "";
    }
  }
  //number of events:
  try {
    LOG("Main", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pNOTICE) << "Using default number of events";
      gOptNevents = kDefOptNevents;
    }
  }
  gOptNevents = TMath::Max(gOptNevents,1); // generate at least 1

  //run number:
  try {
    LOG("Main", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }
  //spline building option:
  gOptBuildSplines =
                genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  //beam direction:
  string direction = "";
  try {
    LOG("Main", pINFO) << "Reading neutrino beam direction";
    direction = genie::utils::clap::CmdLineArgAsString(argc,argv,'d');

    // split the comma separated list
    vector<string> dirv = utils::str::Split(direction, ",");
    assert(dirv.size() == 3);
    double dx = atof( dirv[0].c_str() );
    double dy = atof( dirv[1].c_str() );
    double dz = atof( dirv[2].c_str() );
    // can not pass negative numbers in command line argument as the -
    // symbol marks option names - pass negative number by a trailing m
    dx *= ( (dirv[0].find("m") != string::npos) ? -1. : 1);
    dy *= ( (dirv[1].find("m") != string::npos) ? -1. : 1);
    dz *= ( (dirv[2].find("m") != string::npos) ? -1. : 1);
    gOptBeamDirection.SetXYZ(dx,dy,dz);
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "No input direction - Using default";
      gOptBeamDirection = kDefOptBeamDirection;
    }
  }

  //beam spot:
  string bspot = "";
  try {
    LOG("Main", pINFO) << "Reading neutrino beam spot";
    bspot = genie::utils::clap::CmdLineArgAsString(argc,argv,'b');

    // split the comma separated list
    vector<string> bsv = utils::str::Split(bspot, ",");
    assert(bsv.size() == 3);
    double x = atof( bsv[0].c_str() );
    double y = atof( bsv[1].c_str() );
    double z = atof( bsv[2].c_str() );
    x *= ( (bsv[0].find("m") != string::npos) ? -1. : 1);
    y *= ( (bsv[1].find("m") != string::npos) ? -1. : 1);
    z *= ( (bsv[2].find("m") != string::npos) ? -1. : 1);
    gOptBeamSpot.SetXYZ(x,y,z);
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "No input beam spot - Using default";
      gOptBeamSpot = kDefOptBeamSpot;
    }
  }

  //beam size:
  try {
    LOG("Main", pINFO) << "Reading beam size";
    gOptBeamRt = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "No input beam size - Using default";
      gOptBeamRt = kDefOptBeamRt;
    }
  }
  gOptBeamRt = TMath::Abs(gOptBeamRt); // must be positive

  LOG("Main", pINFO) << "Command line options - Summary:";
  LOG("Main", pINFO) << "run number:           " << gOptRunNu;
  LOG("Main", pINFO) << "spline building:      " << gOptBuildSplines;
  LOG("Main", pINFO) << "number of events:     " << gOptNevents;
  LOG("Main", pINFO) << "detector geom. file:  " << gOptRootGeom;
  LOG("Main", pINFO) << "detector geom. units: " << gOptGeomUnits;
  LOG("Main", pINFO) << "external max pl-list: " << gOptExtMaxPlXml;
  LOG("Main", pINFO) << "beam size:            " << gOptBeamRt;
  LOG("Main", pINFO) << "beam direction:       "
                     << utils::print::Vec3AsString(&gOptBeamDirection);
  LOG("Main", pINFO) << "beam spot:            "
                     << utils::print::Vec3AsString(&gOptBeamSpot);
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("Main", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   testMCJobDriver -f filename [-u units] [-r run] [-n nevents]\n"
    << "                   [-s] [-m filename] [-d dx,dy,dz] [-b x,y,z] [-t size]";
}
//____________________________________________________________________________


