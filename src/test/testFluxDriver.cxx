//____________________________________________________________________________
/*!

\program testFluxDriver

\brief   test program for flux drivers

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 22, 2005
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TF1.h>

#include "FluxDrivers/GCylindTH1Flux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::flux;

const unsigned int kNEvents = 10000;

void testGCylindTH1Flux(void);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("Main", pINFO)  << "Running [GCylindTH1Flux] driver test";

  testGCylindTH1Flux();

  LOG("Main", pINFO)  << "Done!";

  return 0;
}
//___________________________________________________________________
void testGCylindTH1Flux(void)
{
  LOG("Main", pINFO)  << "Creating [GCylindTH1Flux] flux driver";

  GCylindTH1Flux * flux = new GCylindTH1Flux;

  LOG("Main", pINFO)  << "Setting configuration data";

  TF1 * f1 = new TF1("f1","1./x+2.",0.5,5.0);
  TH1D * spectrum1 = new TH1D("spectrum1","numu E",    20,0.5,5);
  spectrum1->FillRandom("f1",100000);

  TF1 * f2 = new TF1("f2","1./(1.+x)",0.5,5.0);
  TH1D * spectrum2 = new TH1D("spectrum2","numubar E", 20,0.5,5);
  spectrum2->FillRandom("f2",1000);

  TVector3 direction(0,0,1);
  TVector3 beam_spot(0,0,-10);

  double Rtransverse = 0.5;

  LOG("Main", pINFO)  << "Configuring [GCylindTH1Flux] flux driver";

  flux -> SetNuDirection      (direction);
  flux -> SetBeamSpot         (beam_spot);
  flux -> SetTransverseRadius (Rtransverse);
  flux -> AddEnergySpectrum   (kPdgNuMu,    spectrum1);
  flux -> AddEnergySpectrum   (kPdgNuMuBar, spectrum2);

  LOG("Main", pINFO) << "Generating flux neutrinos";

  TNtuple * fluxntp = new TNtuple("fluxntp",
           "[GCylindTH1Flux] flux data", "x:y:z:t:px:py:pz:E:pdgc");

  unsigned int ievent = 0;

  while(ievent++ < kNEvents) {

    flux->GenerateNext();

    int pdgc = flux->PdgCode();
    const TLorentzVector & x4 = flux->Position();
    const TLorentzVector & p4 = flux->Momentum();

    fluxntp->Fill( x4.X(),  x4.Y(),  x4.Z(),  x4.T(),
                   p4.Px(), p4.Py(), p4.Pz(), p4.E(), pdgc);
  }

  TFile f("./GCylindTH1Flux.root","recreate");
  fluxntp->Write();
  f.Close();

  delete f1;
  delete f2;
  delete flux;
  delete fluxntp;
}
//___________________________________________________________________

