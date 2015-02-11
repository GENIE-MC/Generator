//____________________________________________________________________________
/*!

\program gtestFluxAtmo

\brief   Test program for atmospheric flux drivers

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created August 22, 2005

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TF1.h>

#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GBartolAtmoFlux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::flux;

const unsigned int kNEvents = 100000;

TNtuple * runGFlukaAtmo3DFluxDriver (void);
TNtuple * runGBartolAtmoFluxDriver  (void);
TNtuple * createFluxNtuple          (GFluxI * flux);

//____________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  LOG("test", pINFO)  << "Running GFlukaAtmo3DFlux driver test";
  TNtuple * ntfluka = runGFlukaAtmo3DFluxDriver();
  ntfluka->SetTitle("GFlukaAtmo3DFlux driver data");

  LOG("test", pINFO)  << "Running GBartolAtmoFlux driver test";
  TNtuple * ntbartol = runGBartolAtmoFluxDriver();
  ntbartol->SetTitle("GBartolAtmoFlux driver data");

  LOG("test", pINFO) << "Saving flux ntuples";

  TFile f("./genie-flux-drivers.root","recreate");
  ntfluka  -> Write("ntfluka");
  ntbartol -> Write("ntbartol");
  f.Close();

  LOG("test", pINFO)  << "Done!";

  return 0;
}
//____________________________________________________________________________
TNtuple * runGFlukaAtmo3DFluxDriver(void)
{
  string base_dir = 
    (gSystem->Getenv("GFLUX_FLUKA3DATMO") ?
       gSystem->Getenv("GFLUX_FLUKA3DATMO") : ".");

  double Rlongitudinal  = 1000.; //m
  double Rtransverse    = 100.;  //m

  GFlukaAtmo3DFlux * flux = new GFlukaAtmo3DFlux;

  LOG("test", pINFO) << base_dir + "/sdave_numu07.dat";
  LOG("test", pINFO) << base_dir + "/sdave_anumu07.dat";
  LOG("test", pINFO) << base_dir + "/sdave_nue07.dat";
  LOG("test", pINFO) << base_dir + "/sdave_anue07.dat";


  flux -> SetFluxFile ( kPdgNuMu,     base_dir + "/sdave_numu07.dat"  );
  flux -> SetFluxFile ( kPdgAntiNuMu, base_dir + "/sdave_anumu07.dat" );
  flux -> SetFluxFile ( kPdgNuE,      base_dir + "/sdave_nue07.dat"   );
  flux -> SetFluxFile ( kPdgAntiNuE,  base_dir + "/sdave_anue07.dat"  );
  flux -> SetRadii(Rlongitudinal, Rtransverse);
  flux -> LoadFluxData();
  flux -> GenerateWeighted(true);
//flux -> ForceMaxEnergy(3);

  LOG("test", pINFO) << "Generating events";
  TNtuple * fluxntp = createFluxNtuple(dynamic_cast<GFluxI*>(flux));
  return fluxntp;
}
//____________________________________________________________________________
TNtuple * runGBartolAtmoFluxDriver(void)
{
  string base_dir = 
    (gSystem->Getenv("GFLUX_BGLRS3DATMO") ?
       gSystem->Getenv("GFLUX_BGLRS3DATMO") : ".");

  double Rlongitudinal  = 1000.; //m
  double Rtransverse    = 100.;  //m

  GBartolAtmoFlux * flux = new GBartolAtmoFlux;

  flux -> SetFluxFile ( kPdgNuMu,     base_dir + "/f210_3_z.kam_num" );
  flux -> SetFluxFile ( kPdgAntiNuMu, base_dir + "/f210_3_z.kam_nbm" );
  flux -> SetFluxFile ( kPdgNuE,      base_dir + "/f210_3_z.kam_nue" );
  flux -> SetFluxFile ( kPdgAntiNuE,  base_dir + "/f210_3_z.kam_nbe" );
  flux -> SetRadii(Rlongitudinal, Rtransverse);
  flux -> LoadFluxData();
  flux -> GenerateWeighted(true);
//flux -> ForceMaxEnergy(300);

  LOG("test", pINFO) << "Generating events";
  TNtuple * fluxntp = createFluxNtuple(dynamic_cast<GFluxI*>(flux));
  delete flux;

  return fluxntp;
}
//____________________________________________________________________________
TNtuple * createFluxNtuple(GFluxI * flux)
{
  TNtuple * fluxntp = 
     new TNtuple("fluxntp", "flux", "x:y:z:t:px:py:pz:E:pdgc:wght");

  LOG("test", pINFO) << "Generating flux neutrinos";

  unsigned int ievent = 0;
  while(ievent++ < kNEvents) {
    LOG("test", pINFO)  << "Event number: " << ievent;
    flux->GenerateNext();
    int    pdgc = flux->PdgCode();
    double wght = flux->Weight();
    const TLorentzVector & x4 = flux->Position();
    const TLorentzVector & p4 = flux->Momentum();
    fluxntp->Fill( 
      x4.X(),  x4.Y(),  x4.Z(),  x4.T(),
      p4.Px(), p4.Py(), p4.Pz(), p4.E(), pdgc, wght);
  }
  return fluxntp;
}
//____________________________________________________________________________

