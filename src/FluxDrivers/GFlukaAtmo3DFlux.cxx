//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 03, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   In vrs 2.3.1, most of the driver implementation code was factored out 
   to the new GAtmoFlux base class so as to be shared by the functionally 
   identical GBartolAtmoFlux driver
 @ Feb 05, 2008 - Chris Backhouse (Oxford)
   Fixed a bug in bin definitions (the TH2 constructor takes an array of 
   bin lower edges, but the bin centres were being passed in instead).
   Added an ad-hoc scaling factor to get agreement with the Bartol flux 
   normalization (needed when both fluxes are stitched together in a single
   simulation with the Bartol flux used as a high energy extension).
*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>

#include <TH2D.h>
#include <TMath.h>

#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "Messenger/Messenger.h"

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GFlukaAtmo3DFlux::GFlukaAtmo3DFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
       << "Instantiating the Fluka-3D atmospheric neutrino flux driver";

  this->SetCosThetaBins(kNGFlk3DCos, kGFlk3DCos);
  this->SetEnergyBins  (kNGFlk3DEv,  kGFlk3DEv);

  this->Initialize();
}
//___________________________________________________________________________
GFlukaAtmo3DFlux::~GFlukaAtmo3DFlux()
{

}
//___________________________________________________________________________
bool GFlukaAtmo3DFlux::FillFluxHisto2D(TH2D * histo, string filename)
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
  char   j1, j2;

  const int kNLines = kNGFlk3DCos * kNGFlk3DEv;
  int iline = 0;
  while (++iline<=kNLines) {
    flux_stream >> energy >> j1 >> costheta >> j2 >> flux;
    LOG("Flux", pDEBUG)
      << "Flux[Ev = " << energy << ", cos8 = " << costheta << "] = " << flux;
    histo->Fill( (Axis_t)energy, (Axis_t)costheta, (Stat_t)flux);
  }
  return true;
}
//___________________________________________________________________________
