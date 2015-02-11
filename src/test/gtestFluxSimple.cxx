//____________________________________________________________________________
/*!

\program gtestFluxSimple

\brief   Test program for simple flux drivers

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

#include "FluxDrivers/GCylindTH1Flux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::flux;

const unsigned int kNEvents = 10000;

TNtuple * runGCylindTH1FluxDriver   (void);
TNtuple * createFluxNtuple          (GFluxI * flux);

//___________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  LOG("test", pINFO)  << "Running GCylindTH1Flux driver test";
  TNtuple * ntcylh1f = runGCylindTH1FluxDriver();
  ntcylh1f->SetTitle("GCylindTH1Flux driver data");

  LOG("test", pINFO) << "Saving flux ntuples";

  TFile f("./genie-flux-drivers.root","recreate");
  ntcylh1f  -> Write("ntcylh1f");
  f.Close();

  delete ntcylh1f;

  LOG("test", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
TNtuple * runGCylindTH1FluxDriver(void)
{
  LOG("test", pINFO)  << "Creating GCylindTH1Flux flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  LOG("test", pINFO)  << "Setting configuration data";

  TF1 * f1 = new TF1("f1","1./x",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu E",    20,0.5,5);
  spectrum1->FillRandom("f1",100000);

  TF1 * f2 = new TF1("f2","x",0.5,5.0);
  TH1D * spectrum2 = new TH1D("spectrum2","numubar E", 20,0.5,5);
  spectrum2->FillRandom("f2",10000);

  TVector3 direction(0,0,1);
  TVector3 beam_spot(0,0,-10);

  double Rtransverse = 0.5;

  LOG("test", pINFO)  << "Configuring GCylindTH1Flux flux driver";

  flux -> SetNuDirection      (direction);
  flux -> SetBeamSpot         (beam_spot);
  flux -> SetTransverseRadius (Rtransverse);
  flux -> AddEnergySpectrum   (kPdgNuMu,     spectrum1);
  flux -> AddEnergySpectrum   (kPdgAntiNuMu, spectrum2);

  LOG("test", pINFO) << "Creating flux ntuple";
  GFluxI * fluxi = dynamic_cast<GFluxI*>(flux);

  TNtuple * fluxntp = createFluxNtuple(fluxi);

  delete f1;
  delete f2;
  delete flux;

  return fluxntp;
}
//___________________________________________________________________
TNtuple * createFluxNtuple(GFluxI * flux)
{
  TNtuple * fluxntp =
      new TNtuple("fluxntp",
            "flux data", "x:y:z:t:px:py:pz:E:pdgc");

  LOG("test", pINFO) << "Generating flux neutrinos";

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

