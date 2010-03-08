//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University

 For the class documentation see the corresponding header file.
 @ Feb 05, 2008 - CB
   This concrete flux driver was added in 2.3.1 by C.Backhouse (Oxford U.)
 @ Feb 23, 2010 - CA
   Build bin arrays at ctor. Re-structuring and clean-up.

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
       << "Instantiating the BGLRS atmospheric neutrino flux driver";

  this->SetBinSizes();
  this->Initialize();
}
//___________________________________________________________________________
GBartolAtmoFlux::~GBartolAtmoFlux()
{

}
//___________________________________________________________________________
void GBartolAtmoFlux::SetBinSizes(void)
{
// Generate the correct cos(theta) and energy bin sizes.
// The flux is given in 20 bins of cos(zenith angle) from -1.0 to 1.0
// (bin width = 0.1) and 30 equally log-spaced energy bins (10 bins
// per decade), with Emin = 10.00 GeV.
//
     
  fCosThetaBins  = new double [kBGLRS3DNumCosThetaBins  + 1];
  fEnergyBins    = new double [kBGLRS3DNumLogEvBins     + 1];
   
  double dcostheta =
      (kBGLRS3DCosThetaMax - kBGLRS3DCosThetaMin) / 
      (double) kBGLRS3DNumCosThetaBins;
     
  double logEmax = TMath::Log10(100.);
  double logEmin = TMath::Log10(kBGLRS3DEvMin);
  double dlogE =
      (logEmax - logEmin) /
      (double) kBGLRS3DNumLogEvBinsPerDecade;
  
  for(unsigned int i=0; i<= kBGLRS3DNumCosThetaBins; i++) {
     fCosThetaBins[i] = kBGLRS3DCosThetaMin + i * dcostheta;
     if(i != kBGLRS3DNumCosThetaBins) {
       LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: CosTheta bin " << i+1
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: CosTheta bin " << kBGLRS3DNumCosThetaBins
         << ": upper edge = " << fCosThetaBins[kBGLRS3DNumCosThetaBins];
     }
  }
      
  for(unsigned int i=0; i<= kBGLRS3DNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kBGLRS3DNumLogEvBins) {
       LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": lower edge = " << fEnergyBins[i];
     } else {
       LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << kBGLRS3DNumLogEvBins
         << ": upper edge = " << fEnergyBins[kBGLRS3DNumLogEvBins];
     }   
  }      

  fNumCosThetaBins = kBGLRS3DNumCosThetaBins;
  fNumEnergyBins   = kBGLRS3DNumLogEvBins; 
}
//___________________________________________________________________________
bool GBartolAtmoFlux::FillFluxHisto2D(TH2D * histo, string filename)
{
  LOG("Flux", pNOTICE) << "Loading: " << filename;

  if(!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
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

  const int kNLines = kBGLRS3DNumCosThetaBins * kBGLRS3DNumLogEvBins;
  int iline = 0;
  while (++iline<=kNLines) {
    flux_stream >> energy >> costheta >> flux >> junkd >> junkd;
    // Compensate for logarithmic units - dlogE=dE/E
    flux /= energy;
    LOG("Flux", pINFO)
      << "Flux[Ev = " << energy 
      << ", cos8 = " << costheta << "] = " << flux;
    histo->Fill( (Axis_t)energy, (Axis_t)costheta, (Stat_t)flux);
  }
  return true;
}
//___________________________________________________________________________
