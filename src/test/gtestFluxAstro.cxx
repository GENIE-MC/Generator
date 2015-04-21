//____________________________________________________________________________
/*!

\program gtestFluxAstro

\brief   Test program for astrophysical neutrino flux drivers

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created March 26, 2010

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "FluxDrivers/GAstroFlux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::flux;

const unsigned int kNEvents = 1000000;

//____________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  const double pi = constants::kPi;

  const double latitude  = pi/5;
  const double longitude = pi/4;
  const double depth     = 3.0*units::km;
  const double size      = 2.0*units::km;

  GDiffuseAstroFlux * difflx = new GDiffuseAstroFlux;

  difflx->ForceMinEnergy(1E+2);
  difflx->ForceMaxEnergy(1E+8);
  difflx->GenerateWeighted(true);
  difflx->SetDetectorPosition(latitude,longitude,depth,size);
  difflx->SetRelNuPopulations(1,2,0,1,2,0);
  difflx->SetEnergyPowLawIdx(3.5);

  TNtuple * fluxntp = 
     new TNtuple("fluxntp", "flux", "x:y:z:t:px:py:pz:E:pdgc:wght");

  unsigned int ievent = 0;
  while(ievent++ < kNEvents) {
    LOG("test", pINFO)  << "Event number: " << ievent;
    difflx->GenerateNext();
    int    pdgc = difflx->PdgCode();
    double wght = difflx->Weight();
    const TLorentzVector & x4 = difflx->Position();
    const TLorentzVector & p4 = difflx->Momentum();
    fluxntp->Fill( 
      x4.X(),  x4.Y(),  x4.Z(),  x4.T(),
      p4.Px(), p4.Py(), p4.Pz(), p4.E(), pdgc, wght);
  }

  TFile f("./genie-astro-flux.root","recreate");
  fluxntp->Write();
  f.Close();

  LOG("test", pINFO)  << "Done!";

  return 0;
}
//____________________________________________________________________________
