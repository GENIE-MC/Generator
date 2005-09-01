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

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJob.h"
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using std::string;
using namespace genie;
using namespace genie::flux;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- Default geometry
  string base_dir = string( gSystem->Getenv("GENIE") );
  string filename = base_dir+ string("/src/test/TestGeometry.root");
  //-- Scan for filename from the command line argument (following -f)
  for(int iarg = 0; iarg < argc-1; iarg++) {
     string argument(argv[iarg]);
     if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }

  //-- Create the GENIE MC-job driver
  GMCJob mcj;

  //-- Specify a flux driver
/*
  LOG("Main", pINFO)  << "Creating [GFlukaAtmo3DFlux] flux driver";

  GFlukaAtmo3DFlux * flux = new GFlukaAtmo3DFlux;

  flux->SetNuMuFluxFile("/home/costas/tmp/downloads/sdave_numu07.dat");
  flux->SetNuMuBarFluxFile("/home/costas/tmp/downloads/sdave_anumu07.dat");
  flux->SetNuEFluxFile("/home/costas/tmp/downloads/sdave_nue07.dat");
  flux->SetNuEBarFluxFile("/home/costas/tmp/downloads/sdave_anue07.dat");
  flux->SetRadii(1000.,100.);

  flux->LoadFluxData();
*/

  LOG("Main", pINFO)  << "Creating [GCylindTH1Flux] flux driver";

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

  //-- Specify the geometry analyzer

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(filename);

  //-- Set the flux and the geometry analyzer to the GENIE MC driver

  LOG("Main", pINFO)
    << "Creating the GENIE MC Job Driver & specifying flux & geometry";

  GFluxI *        fluxb = dynamic_cast<GFluxI *>       (flux);
  GeomAnalyzerI * geomb = dynamic_cast<GeomAnalyzerI *>(geom);

  mcj.UseFluxDriver  (fluxb);
  mcj.UseGeomAnalyzer(geomb);

  //-- Configure the GENIE MC driver

  mcj.Configure();

  //-- Start generating events -here, just 1 for testing purposes-

  EventRecord * event = mcj.GenerateEvent();

  LOG("Main", pINFO) << *event;

  delete f1;
  delete flux;
  delete geom;

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
