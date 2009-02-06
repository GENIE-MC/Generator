//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University, January 26, 2008

 For the class documentation see the corresponding header file.
 @ Feb 05, 2008 - CB
   This concrete flux driver was added in 2.3.1 by C.Backhouse (Oxford U.)
*/
//____________________________________________________________________________

#include <fstream>
#include <cassert>

#include <TH2D.h>
#include <TMath.h>

#include "FluxDrivers/GBartolAtmoFlux.h"
#include "Messenger/Messenger.h"

using std::ifstream;
using std::ios;

using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GBartolAtmoFlux::GBartolAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
       << "Instantiating the Bartol atmospheric neutrino flux driver";

  this->SetCosThetaBins(kNGBrtCos, kGBrtCos);
  this->SetEnergyBins  (kNGBrtEv,  kGBrtEv);

  this->Initialize();
}
//___________________________________________________________________________
GBartolAtmoFlux::~GBartolAtmoFlux()
{

}
//___________________________________________________________________________
bool GBartolAtmoFlux::FillFluxHisto2D(TH2D * histo, string filename)
{
  LOG("Flux", pNOTICE) << "Loading: " << filename;

  if(!histo) {
     LOG("Flux", pERROR) << "The flux histo to fill is not instantiated!";
     return false;
  }
  if(filename.size() == 0) {
     // the user wants to skip this flux component - create an empty flux
     this->ZeroFluxHisto2D(histo);
     fNSkipped++;
     assert(fNSkipped<4); // you must load at least one flux file
     return true;
  }

  ifstream flux_stream(filename.c_str(), ios::in);
  if(!flux_stream) {
     LOG("Flux", pERROR) << "Could not open file: " << filename;
     return false;
  }

  double energy, costheta, flux;
  double junkd; // throw away error estimates

  // throw away comment line
  flux_stream.ignore(99999, '\n');

  const int kNLines = kNGBrtCos * kNGBrtEv;
  int iline = 0;
  while (++iline<=kNLines) {
    flux_stream >> energy >> costheta >> flux >> junkd >> junkd;
    // Compensate for logarithmic units - dlogE=dE/E
    flux /= energy;
    LOG("Flux", pDEBUG)
      << "Flux[Ev = " << energy << ", cos8 = " << costheta << "] = " << flux;
    histo->Fill( (Axis_t)energy, (Axis_t)costheta, (Stat_t)flux);
  }
  return true;
}
//___________________________________________________________________________
