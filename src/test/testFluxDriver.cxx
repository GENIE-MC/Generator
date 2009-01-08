//____________________________________________________________________________
/*!

\program gtestFluxDriver

\brief   Test program for flux drivers

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 22, 2005

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TF1.h>

#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::flux;

const unsigned int kNEvents = 10000;

TNtuple * runGCylindTH1FluxDriver   (void);
TNtuple * runGFlukaAtmo3DFluxDriver (void);
TNtuple * createFluxNtuple           (GFluxI * flux);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("Main", pINFO)  << "Running [GCylindTH1Flux] driver test";
  TNtuple * ntcylf = runGCylindTH1FluxDriver();
  ntcylf->SetTitle("[GCylindTH1Flux] driver data");

  LOG("Main", pINFO)  << "Running [GFlukaAtmo3DFlux] driver test";
  TNtuple * ntflk3d = runGFlukaAtmo3DFluxDriver();
  ntflk3d->SetTitle("[GFlukaAtmo3DFlux] driver data");

  LOG("Main", pINFO) << "Saving flux ntuples";

  TFile f("./genie-flux-drivers.root","recreate");
  ntcylf  -> Write("ntcylf");
  ntflk3d -> Write("ntflk3d");
  f.Close();

  delete ntcylf;
  delete ntflk3d;

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
TNtuple * runGCylindTH1FluxDriver(void)
{
  LOG("Main", pINFO)  << "Creating [GCylindTH1Flux] flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  LOG("Main", pINFO)  << "Setting configuration data";

  TF1 * f1 = new TF1("f1","1./x",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu E",    20,0.5,5);
  spectrum1->FillRandom("f1",100000);

  TF1 * f2 = new TF1("f2","x",0.5,5.0);
  TH1D * spectrum2 = new TH1D("spectrum2","numubar E", 20,0.5,5);
  spectrum2->FillRandom("f2",10000);

  TVector3 direction(0,0,1);
  TVector3 beam_spot(0,0,-10);

  double Rtransverse = 0.5;

  LOG("Main", pINFO)  << "Configuring [GCylindTH1Flux] flux driver";

  flux -> SetNuDirection      (direction);
  flux -> SetBeamSpot         (beam_spot);
  flux -> SetTransverseRadius (Rtransverse);
  flux -> AddEnergySpectrum   (kPdgNuMu,     spectrum1);
  flux -> AddEnergySpectrum   (kPdgAntiNuMu, spectrum2);

  LOG("Main", pINFO) << "Creating flux ntuple";
  GFluxI * fluxi = dynamic_cast<GFluxI*>(flux);

  TNtuple * fluxntp = createFluxNtuple(fluxi);

  delete f1;
  delete f2;
  delete flux;

  return fluxntp;
}
//___________________________________________________________________
TNtuple * runGFlukaAtmo3DFluxDriver(void)
{
  LOG("Main", pINFO)  << "Creating [GFlukaAtmo3DFlux] flux driver";

  GFlukaAtmo3DFlux * flux = new GFlukaAtmo3DFlux;

  LOG("Main", pINFO)  << "Setting configuration data";

  string base_dir = (gSystem->Getenv("GFLUX_FLUKA3DATMO") ?
                        gSystem->Getenv("GFLUX_FLUKA3DATMO") : ".");

  string numu_flux_file    = base_dir + "/sdave_numu07.dat";
  string numubar_flux_file = base_dir + "/sdave_anumu07.dat";
  string nue_flux_file     = base_dir + "/sdave_nue07.dat";
  string nuebar_flux_file  = base_dir + "/sdave_anue07.dat";
  double Rlongitudinal     = 1000.; //m
  double Rtransverse       = 100.;  //m

  LOG("Main", pINFO)  << "Configuring [GFlukaAtmo3DFlux] flux driver";

  flux -> SetNuMuFluxFile    ( numu_flux_file    );
  flux -> SetNuMuBarFluxFile ( numubar_flux_file );
  flux -> SetNuEFluxFile     ( nue_flux_file     );
  flux -> SetNuEBarFluxFile  ( nuebar_flux_file  );

  flux -> SetRadii(Rlongitudinal, Rtransverse);
  flux -> LoadFluxData();

  LOG("Main", pINFO) << "Creating flux ntuple";

  GFluxI * fluxi = dynamic_cast<GFluxI*>(flux);

  TNtuple * fluxntp = createFluxNtuple(fluxi);

  delete flux;

  return fluxntp;
}
//___________________________________________________________________
TNtuple * createFluxNtuple(GFluxI * flux)
{
  TNtuple * fluxntp = new TNtuple("fluxntp",
                             "flux data", "x:y:z:t:px:py:pz:E:pdgc");

  LOG("Main", pINFO) << "Generating flux neutrinos";

  unsigned int ievent = 0;
  while(ievent++ < kNEvents) {

    flux->GenerateNext();

    int pdgc = flux->PdgCode();
    const TLorentzVector & x4 = flux->Position();
    const TLorentzVector & p4 = flux->Momentum();

    fluxntp->Fill( x4.X(),  x4.Y(),  x4.Z(),  x4.T(),
                   p4.Px(), p4.Py(), p4.Pz(), p4.E(), pdgc);
  }
  return fluxntp;
}
//___________________________________________________________________

