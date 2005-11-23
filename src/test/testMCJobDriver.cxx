//____________________________________________________________________________
/*!

\program testMCJobDriver

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 22, 2005
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "PDG/PDGCodes.h"
#include "Utils/XSecSplineList.h"
#include "Utils/UnitUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

void GetCommandLineArgs(int argc, char ** argv);

//command line options
bool   gOptBuildSplines; // spline building option
int    gOptNevents;      // number of events to generate
string gOptRootGeom;     // detector geometry ROOT file
string gOptGeomUnits;    // detector geometry units

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- Parse command line arguments
  GetCommandLineArgs(argc, argv);

  //-- Create/configure a flux driver
  LOG("Main", pINFO)
            << "Creating/configuring the GCylindTH1Flux flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  TF1 * f1 = new TF1("f1","1./x+2.",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu spectrum", 20,0.5,5);
  spectrum1->FillRandom("f1",10000);

  TVector3 direction(0,0,1);
  TVector3 beam_spot(0,0,-10);

  flux -> SetNuDirection      (direction);
  flux -> SetBeamSpot         (beam_spot);
  flux -> SetTransverseRadius (0.5);
  flux -> AddEnergySpectrum   (kPdgNuMu, spectrum1);

  //-- Create/configure a geometry driver
  LOG("Main", pINFO) << "Creating/configuring the ROOT geom. driver";

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptRootGeom);
  geom->SetUnits(genie::utils::units::UnitFromString(gOptGeomUnits));

  //-- Create the GENIE MC-job driver
  LOG("Main", pINFO) << "Creating the GENIE MC Job driver";

  GMCJDriver mcj;

  //-- Load the the flux and the geometry drivers to the MC driver

  LOG("Main", pINFO)
         << "Loading flux & geometry drivers to GENIE MC Job driver";

  GFluxI *        fluxb = dynamic_cast<GFluxI *>       (flux);
  GeomAnalyzerI * geomb = dynamic_cast<GeomAnalyzerI *>(geom);

  mcj.UseFluxDriver  (fluxb);
  mcj.UseGeomAnalyzer(geomb);

  //-- Configure the GENIE MC driver
  LOG("Main", pINFO) << "Configuring the GENIE MC Job driver";

  mcj.Configure();

  //-- If this job uses cross section splines, build all splines that
  //   are needed and have not already loaded from an XML file via
  //   XSecSplineList::AutoLoad()
  if(gOptBuildSplines) mcj.UseSplines();

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kNFEventRecord);
  ntpw.InitTree("mcjobdriver.root");

  //-- Start generating events

  int i=0;
  while (i<gOptNevents) {
     EventRecord * event = mcj.GenerateEvent();

     LOG("Main", pINFO) << *event;

     ntpw.AddEventRecord(i++, event);
     delete event;
  }

  //-- save the ntuple
  ntpw.SaveTree();

  delete f1;
  delete flux;
  delete geom;

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // default options
  string kDefOptRootGeom  = string(gSystem->Getenv("GENIE")) + 
                            "/src/test/data/GeometryLArPbBox.root";
  string kDefOptGeomUnits = "m";
  int    kDefOptNevents   = 10;

  //geometry file:
  try {
    LOG("Main", pINFO) << "Getting input geometry file";
    gOptRootGeom =
              genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pNOTICE) << "Using default geometry file";
      gOptRootGeom = kDefOptRootGeom;
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

  //spline building option:
  gOptBuildSplines =
                genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  LOG("Main", pINFO) << "Command line options - Summary:";
  LOG("Main", pINFO) << "spline building:      " << gOptBuildSplines;
  LOG("Main", pINFO) << "number of events:     " << gOptNevents;
  LOG("Main", pINFO) << "detector geom. file:  " << gOptRootGeom;
  LOG("Main", pINFO) << "detector geom. units: " << gOptGeomUnits;
}
//___________________________________________________________________

