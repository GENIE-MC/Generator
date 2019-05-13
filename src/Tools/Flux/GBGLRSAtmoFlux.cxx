//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Christopher Backhouse <c.backhouse1@physics.ox.ac.uk>
         Oxford University

 For the class documentation see the corresponding header file.
 @ Feb 05, 2008 - CB
   This concrete flux driver was added in 2.3.1 by C.Backhouse (Oxford U.)
 @ Feb 23, 2010 - CA
   Build bin arrays at ctor. Re-structuring and clean-up.
 @ Feb 23, 2012 - AB
   Combine the flux calculations at low and high energies.

*/
//____________________________________________________________________________

#include <fstream>
#include <cassert>

#include <TH3D.h>
#include <TMath.h>

#include "Tools/Flux/GBGLRSAtmoFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GBGLRSAtmoFlux,genie::flux::GBGLRSAtmoFlux)

using std::ifstream;
using std::ios;

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
GBGLRSAtmoFlux::GBGLRSAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
     << "Instantiating the GENIE BGLRS atmospheric neutrino flux driver";

  this->Initialize();
  this->SetBinSizes();
}
//___________________________________________________________________________
GBGLRSAtmoFlux::~GBGLRSAtmoFlux()
{

}
//___________________________________________________________________________
void GBGLRSAtmoFlux::SetBinSizes(void)
{
// Generate the correct cos(theta) and energy bin sizes.
//
// Zenith angle binning: the flux is given in 20 bins of 
// cos(zenith angle) from -1.0 to 1.0 (bin width = 0.1) 
//
// Neutrino energy binning: the Bartol flux files are 
// provided in two pieces 
//  (1) low energy piece (<10 GeV), solar min or max,
//      given in 40 log-spaced bins from 0.1 to 10 GeV 
//      (20 bins per decade)
//  (2) high energy piece (>10 GeV), without solar effects, 
//      given in 30 log-spaced bins from 10 to 1000 GeV
//      (10 bins per decade)
     
  fPhiBins       = new double [2];
  fCosThetaBins  = new double [kBGLRS3DNumCosThetaBins + 1];
  fEnergyBins    = new double [kBGLRS3DNumLogEvBinsLow + kBGLRS3DNumLogEvBinsHigh + 1];

  fPhiBins[0] = 0;
  fPhiBins[1] = 2.*kPi;
   
  double dcostheta =
      (kBGLRS3DCosThetaMax - kBGLRS3DCosThetaMin) / 
      (double) kBGLRS3DNumCosThetaBins;
     
  double logEmin = TMath::Log10(kBGLRS3DEvMin);
  double dlogElow = 1.0 / (double) kBGLRS3DNumLogEvBinsPerDecadeLow;
  double dlogEhigh = 1.0 / (double) kBGLRS3DNumLogEvBinsPerDecadeHigh;

  double costheta = kBGLRS3DCosThetaMin;

  for(unsigned int i=0; i<= kBGLRS3DNumCosThetaBins; i++) {
    if( i==0 ) ; // do nothing
    else costheta += dcostheta;
    fCosThetaBins[i] = costheta;
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
     
  double logE = logEmin;
 
  for(unsigned int i=0; i<=kBGLRS3DNumLogEvBinsLow+kBGLRS3DNumLogEvBinsHigh; i++) {
    if( i==0 ) ; // do nothing
    else if( i<=kBGLRS3DNumLogEvBinsLow ) logE += dlogElow;
    else                                  logE += dlogEhigh;
    fEnergyBins[i] = TMath::Power(10.0, logE);
    if(i != kBGLRS3DNumLogEvBinsLow+kBGLRS3DNumLogEvBinsHigh) {
      LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": lower edge = " << fEnergyBins[i];
    } else {
      LOG("Flux", pDEBUG)
         << "FLUKA 3d flux: Energy bin " << kBGLRS3DNumLogEvBinsLow+kBGLRS3DNumLogEvBinsHigh
         << ": upper edge = " << fEnergyBins[kBGLRS3DNumLogEvBinsLow+kBGLRS3DNumLogEvBinsHigh];
    } 
  }      

  fNumPhiBins      = 1;
  fNumCosThetaBins = kBGLRS3DNumCosThetaBins;
  fNumEnergyBins   = kBGLRS3DNumLogEvBinsLow + kBGLRS3DNumLogEvBinsHigh; 
  fMaxEv = fEnergyBins[fNumEnergyBins];
}
//___________________________________________________________________________
bool GBGLRSAtmoFlux::FillFluxHisto(int nu_pdg, string filename)
{
  LOG("Flux", pNOTICE) 
    << "Loading BGLRS flux for neutrino: " << nu_pdg 
    << " from file: " << filename;

  TH3D* histo = 0;
  std::map<int,TH3D*>::iterator myMapEntry = fRawFluxHistoMap.find(nu_pdg);
  if( myMapEntry != fRawFluxHistoMap.end() ){
      histo = myMapEntry->second;
  }
  if(!histo) {
     LOG("Flux", pERROR) << "Null flux histogram!";
     return false;
  }
  ifstream flux_stream(filename.c_str(), ios::in);
  if(!flux_stream) {
     LOG("Flux", pERROR) << "Could not open file: " << filename;
     return false;
  }

  int ibin;
  double energy, costheta, flux;
  double junkd; // throw away error estimates

  double scale = 1.0; // 1.0 [m^2], OR 1.0e-4 [cm^2]

  // throw away comment line
  flux_stream.ignore(99999, '\n');

  while ( !flux_stream.eof() ) {
    flux = 0.0;
    flux_stream >> energy >> costheta >> flux >> junkd >> junkd;
    if( flux>0.0 ){
      // Compensate for logarithmic units - dlogE=dE/E
      // [Note: should do this explicitly using bin widths]
      flux /= energy;
      LOG("Flux", pINFO)
        << "Flux[Ev = " << energy 
        << ", cos8 = " << costheta << "] = " << flux;
      ibin = histo->FindBin( (Axis_t)energy, (Axis_t)costheta, (Axis_t)kPi );
      histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
    }
  }
  return true;
}
//___________________________________________________________________________
