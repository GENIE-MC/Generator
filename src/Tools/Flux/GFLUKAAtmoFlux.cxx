//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   In vrs 2.3.1, most of the driver implementation code was factored out 
   to the new GAtmoFlux base class so as to be shared by the functionally 
   identical GBGLRSAtmoFlux driver
 @ Feb 05, 2008 - Chris Backhouse (Oxford)
   Fixed a bug in bin definitions (the TH2 constructor takes an array of 
   bin lower edges, but the bin centres were being passed in instead).
 @ Feb 23, 2010 - CA
   Build bin arrays at ctor. Re-structuring and clean-up.

*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>

#include <TH3D.h>
#include <TMath.h>

#include "Tools/Flux/GFLUKAAtmoFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GFLUKAAtmoFlux,genie::flux::GFLUKAAtmoFlux)

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
GFLUKAAtmoFlux::GFLUKAAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
    << "Instantiating the GENIE FLUKA atmospheric neutrino flux driver";

  this->Initialize();
  this->SetBinSizes();
}
//___________________________________________________________________________
GFLUKAAtmoFlux::~GFLUKAAtmoFlux()
{

}
//___________________________________________________________________________
void GFLUKAAtmoFlux::SetBinSizes(void)
{
// Generate the correct cos(theta) and energy bin sizes
// The flux is given in 40 bins of cos(zenith angle) from -1.0 to 1.0
// (bin width = 0.05) and 61 equally log-spaced energy bins (20 bins 
// per decade), with Emin = 0.100 GeV.
//

  fPhiBins       = new double [2];
  fCosThetaBins  = new double [kGFlk3DNumCosThetaBins  + 1];
  fEnergyBins    = new double [kGFlk3DNumLogEvBins     + 1];

  fPhiBins[0] = 0;
  fPhiBins[1] = 2.*kPi;

  double dcostheta = 
      (kGFlk3DCosThetaMax - kGFlk3DCosThetaMin) /
      (double) kGFlk3DNumCosThetaBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(kGFlk3DEvMin);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGFlk3DNumLogEvBinsPerDecade;

  for(unsigned int i=0; i<= kGFlk3DNumCosThetaBins; i++) {
     fCosThetaBins[i] = kGFlk3DCosThetaMin + i * dcostheta;
     if(i != kGFlk3DNumCosThetaBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: CosTheta bin " << kGFlk3DNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGFlk3DNumCosThetaBins];
     }
  }

  for(unsigned int i=0; i<= kGFlk3DNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kGFlk3DNumLogEvBins) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1 
         << ": lower edge = " << fEnergyBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << kGFlk3DNumLogEvBins 
         << ": upper edge = " << fEnergyBins[kGFlk3DNumLogEvBins];
     }
  }

  for(unsigned int i=0; i< kGFlk3DNumLogEvBins; i++) {
       LOG("Flux", pDEBUG) 
         << "FLUKA 3d flux: Energy bin " << i+1
         << ": bin centre = " << (fEnergyBins[i] + fEnergyBins[i+1])/2.;
  }

  fNumPhiBins      = 1;
  fNumCosThetaBins = kGFlk3DNumCosThetaBins;
  fNumEnergyBins   = kGFlk3DNumLogEvBins;
  fMaxEv           = fEnergyBins[fNumEnergyBins];
}
//____________________________________________________________________________
bool GFLUKAAtmoFlux::FillFluxHisto(int nu_pdg, string filename)
{
  LOG("Flux", pNOTICE)
    << "Loading FLUKA atmospheric flux for neutrino: " << nu_pdg 
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

  int    ibin;
  double energy, costheta, flux;
  char   j1, j2;

  double scale = 1.0; // 1.0 [m^2], OR 1.0e-4 [cm^2]

  while ( !flux_stream.eof() ) {
    flux = 0.0;
    flux_stream >> energy >> j1 >> costheta >> j2 >> flux;
    if( flux>0.0 ){
      LOG("Flux", pINFO)
        << "Flux[Ev = " << energy 
        << ", cos8 = " << costheta << "] = " << flux;
      // note: reversing the Fluka sign convention for zenith angle
      //       1 phi bin
      ibin = histo->FindBin( (Axis_t)energy, (Axis_t)(-costheta), (Axis_t)kPi );   
      histo->SetBinContent( ibin, (Stat_t)(scale*flux) );
    }
  }
  return true;
}
//___________________________________________________________________________
