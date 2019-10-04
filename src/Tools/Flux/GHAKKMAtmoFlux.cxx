//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GHAKKMAtmoFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GHAKKMAtmoFlux,genie::flux::GHAKKMAtmoFlux)

using std::ifstream;
using std::getline;
using std::ios;
using namespace genie;
using namespace genie::flux;

//____________________________________________________________________________
GHAKKMAtmoFlux::GHAKKMAtmoFlux() :
GAtmoFlux()
{
  LOG("Flux", pNOTICE)
    << "Instantiating the GENIE HAKKM atmospheric neutrino flux driver";

  this->Initialize();
  this->SetBinSizes();
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
         << "HAKKM 3D flux: CosTheta bin " << i+1 
         << ": lower edge = " << fCosThetaBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "HAKKM 3D flux: CosTheta bin " << kGHnd3DNumCosThetaBins 
         << ": upper edge = " << fCosThetaBins[kGHnd3DNumCosThetaBins];
     }
  }

  //
  // phi
  //

  fPhiBins    = new double [kGHnd3DNumPhiBins + 1];
  fNumPhiBins = kGHnd3DNumPhiBins;

  double d2r = constants::kPi/180.;

  double dphi = 
      d2r * (kGHnd3DPhiMax - kGHnd3DPhiMin) /
      (double) kGHnd3DNumPhiBins;

  for(unsigned int i=0; i<= kGHnd3DNumPhiBins; i++) {
     fPhiBins[i] = kGHnd3DPhiMin + i * dphi;
     if(i != kGHnd3DNumPhiBins) {
       LOG("Flux", pDEBUG) 
         << "HAKKM 3D flux: Phi bin " << i+1 
         << ": lower edge = " << fPhiBins[i];
     } else {
       LOG("Flux", pDEBUG) 
         << "HAKKM 3D flux: Phi bin " << kGHnd3DNumPhiBins 
         << ": upper edge = " << fPhiBins[kGHnd3DNumPhiBins];
     }
  }

  //
  // log(E)
  //

  // For each costheta,phi pair there are N logarithmically spaced
  // neutrino energy values (starting at 0.1 GeV with 20 values per decade
  // up to 10000 GeV) each with corresponding flux values.
  // To construct a flux histogram, use N+1 bins from 0 up to maximum
  // value. Assume that the flux value given for E=E0 is the flux value
  // at the bin that has E0 as its upper edge.
  //
  fEnergyBins    = new double [kGHnd3DNumLogEvBins + 1]; // bin edges
  fNumEnergyBins = kGHnd3DNumLogEvBins;

  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(0.1);
  double dlogE = 
      (logEmax - logEmin) / 
      (double) kGHnd3DNumLogEvBinsPerDecade;

  fEnergyBins[0] = 0;
  for(unsigned int i=0; i < fNumEnergyBins; i++) {
     fEnergyBins[i+1] = TMath::Power(10., logEmin + i*dlogE);
     if(i < kGHnd3DNumLogEvBins) {
       LOG("Flux", pDEBUG) 
          << "HAKKM 3D flux: Energy bin " << i+1 
          << ": upper edge = " << fEnergyBins[i+1] << " GeV";
     } 
  }

  fMaxEv = fEnergyBins[fNumEnergyBins];
  LOG("Flux", pDEBUG) 
    << "HAKKM 3D flux: Maximum energy = " << fMaxEv;

}
//____________________________________________________________________________
bool GHAKKMAtmoFlux::FillFluxHisto(int nu_pdg, string filename)
{
  LOG("Flux", pNOTICE)
    << "Loading HAKKM flux for neutrino: " << nu_pdg
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

  double r2d = 180./constants::kPi;

  int icostheta = kGHnd3DNumCosThetaBins; // starting from last bin [0.9 - 1.0]
  int iphi      = 0;

  while (!flux_stream.eof()) {

    string comment = "";
    getline(flux_stream,comment);
    LOG("Flux", pDEBUG) << "Comment line from HAKKM input file: " << comment;
    getline(flux_stream,comment);
    LOG("Flux", pDEBUG) << "Comment line from HAKKM input file: " << comment;

    iphi++;
    if(iphi == kGHnd3DNumPhiBins+1) {
      icostheta--;
      iphi=1;
    }

    LOG("Flux", pDEBUG)
       << "icostheta = " << icostheta << ", iphi = " << iphi << " / "
       << "costheta bins = " << kGHnd3DNumCosThetaBins << ", phi bins = " << kGHnd3DNumPhiBins;
    LOG("Flux", pDEBUG)
       << "The following set of (energy,flux) values corresponds to "
       << "costheta = [" << fCosThetaBins[icostheta-1] << ", " << fCosThetaBins[icostheta] << "]"
       << ", phi = [" << fPhiBins[iphi-1] << ", " << fPhiBins[iphi] << "] rad"
       << " or [" << r2d * fPhiBins[iphi-1] << ", " << r2d * fPhiBins[iphi] << "] deg) ";

    for(unsigned int i=0; i < kGHnd3DNumLogEvBins; i++) {    

      int ienergy = i+1;

      double energy   = 0;
      double fnumu    = 0;
      double fnumubar = 0;
      double fnue     = 0;
      double fnuebar  = 0;

      flux_stream >> energy >> fnumu >> fnumubar >> fnue >> fnuebar;

      // fitting this easily into what is done for FLUKA, BGLRS where a 
      // different file is specified for each neurtino species means that
      // the input file for HAKKM has to be read 4 times (at most).
      // However, this maintains the ability to switch off individual 
      // components at source and generate interactions for some species only
    
      double flux = 0.;
      if(nu_pdg == kPdgNuMu    ) flux = fnumu;
      if(nu_pdg == kPdgAntiNuMu) flux = fnumubar;
      if(nu_pdg == kPdgNuE     ) flux = fnue;
      if(nu_pdg == kPdgAntiNuE ) flux = fnuebar;
      LOG("Flux", pDEBUG)
         << "Flux (nu_pdg = " << nu_pdg 
         << "; Ev = " << energy << " GeV / bin used = ["
         << fEnergyBins[ienergy-1] << ", " << fEnergyBins[ienergy] << "] GeV"
         << ") = " << flux << " (m^2 sec sr GeV)^-1";
      if(flux > 0.) {
        histo->SetBinContent(ienergy,icostheta,iphi,flux);
      }
    }
    getline(flux_stream,comment);
    LOG("Flux", pDEBUG) << comment;
  }
  return true;
}
//___________________________________________________________________________
