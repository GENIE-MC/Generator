//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <cassert>
#include <fstream>

#include <TH3D.h>
#include <TMath.h>

#include "FluxDrivers/GHAKKMAtmoFlux.h"
#include "Messenger/Messenger.h"

#include "FluxDrivers/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GHAKKMAtmoFlux,genie::flux::GHAKKMAtmoFlux)

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GHAKKMAtmoFlux::GHAKKMAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
    << "Instantiating the ATMNC 3D atmospheric neutrino flux driver";

  this->SetBinSizes();
  this->Initialize();
}
//___________________________________________________________________________
GHAKKMAtmoFlux::~GHAKKMAtmoFlux()
{

}
//___________________________________________________________________________
void GHAKKMAtmoFlux::SetBinSizes(void)
{ 
  //
  // cos(theta)
  //

  fCosThetaBins    = new double [kGHnd3DNumCosThetaBins + 1];
  fNumCosThetaBins = kGHnd3DNumCosThetaBins;

  double dcostheta = 
      (kGHnd3DCosThetaMax - kGHnd3DCosThetaMin) /
      (double) kGHnd3DNumCosThetaBins;

  for(unsigned int i=0; i<= kGHnd3DNumCosThetaBins; i++) {
     fCosThetaBins[i] = kGHnd3DCosThetaMin + i * dcostheta;
     if(i != kGHnd3DNumCosThetaBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: CosTheta bin " << kGHnd3DNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGHnd3DNumCosThetaBins];
     }
  }


  //
  // phi
  //

  fPhiBins    = new double [kGHnd3DNumPhiBins + 1];
  fNumPhiBins = kGHnd3DNumPhiBins;

  double dphi = 
      (kGHnd3DPhiMax - kGHnd3DPhiMin) /
      (double) kGHnd3DNumPhiBins;

  for(unsigned int i=0; i<= kGHnd3DNumPhiBins; i++) {
     fPhiBins[i] = kGHnd3DPhiMin + i * dphi;
     if(i != kGHnd3DNumPhiBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Phi bin " << i+1 
         << ": lower edge = " << fPhiBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Phi bin " << kGHnd3DNumPhiBins 
         << ": upper edge = " << fPhiBins[kGHnd3DNumPhiBins];
     }
  }

  //
  // log(E)
  //

  fEnergyBins    = new double [kGHnd3DNumLogEvBins + 1];
  fNumEnergyBins = kGHnd3DNumLogEvBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(kGHnd3DEvMin);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGHnd3DNumLogEvBinsPerDecade;

  for(unsigned int i=0; i<= kGHnd3DNumLogEvBins; i++) {
     fEnergyBins[i] = TMath::Power(10., logEmin + i*dlogE);
     if(i != kGHnd3DNumLogEvBins) {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Energy bin " << i+1 
         << ": lower edge = " << fEnergyBins[i] << " GeV"
         << ", bin centre = " << (fEnergyBins[i] + fEnergyBins[i+1])/2. << " GeV";
     } else {
       LOG("Flux", pDEBUG) 
         << "ATMNC 3D flux: Energy bin " << kGHnd3DNumLogEvBins 
         << ": upper edge = " << fEnergyBins[kGHnd3DNumLogEvBins] << " GeV";
     }
  }

}
//____________________________________________________________________________
bool GHAKKMAtmoFlux::FillFluxHisto(TH3D * histo, string filename)
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

  double energy, costheta, phi, flux;
  char   j1, j2, j3;

  const int kNLines = kGHnd3DNumCosThetaBins * kGHnd3DNumPhiBins * kGHnd3DNumLogEvBins;
  int iline = 0;
  while (++iline<=kNLines) {
    flux_stream >> energy >> j1 >> costheta >> j2 >> phi >> j3 >> flux;
    LOG("Flux", pINFO)
      << "Flux[Ev = " << energy 
      << ", cos8 = " << costheta 
      << ", phi = " << phi << "] = " << flux;
    histo->Fill( (Axis_t)energy, (Axis_t)costheta, (Axis_t)phi, (Axis_t)flux );
  }
  return true;
}
//___________________________________________________________________________
