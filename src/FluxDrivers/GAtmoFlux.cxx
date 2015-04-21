//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   This class was added in 2.3.1 by code factored out from the concrete
   GFlukaAtmo3DFlux driver & the newer, largely similar, GBartolAtmoFlux
   driver.
 @ Feb 23, 2010 - CA
   Re-structuring and clean-up. Added option to generate weighted flux.
   Added option to specify a maximum energy cut.
 @ Feb 24, 2010 - CA
   Added option to specify a minimum energy cut.
 @ Sep 22, 2010 - TF, CA
   Added SetUserCoordSystem(TRotation &) to specify a rotation from the
   Topocentric Horizontal (THZ) coordinate system to a user-defined 
   topocentric coordinate system. Added NFluxNeutrinos() to get number of
   flux neutrinos generated for sample normalization purposes (note that, in 
   the presence of cuts, this is not the same as the number of flux neutrinos 
   thrown towards the geometry).
 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear and GFluxI::Index as 
   these methods needed for pre-generation of flux interaction probabilities 
   in GMCJDriver.
@ Feb 03, 2011 - TF
   Bug fixed: events are now generated randomly and uniformly on a disc with 
   R = R_{transverse}
@ Feb 23, 2012 - AB
   Bug fixed: events were being generated according to the differential flux
   in each energy bin, dPhi/dE, rather than the total flux, Phi, in each bin.
   This has now been fixed.
*/
//____________________________________________________________________________

#include <cassert>

#include <TH2D.h>
#include <TMath.h>

#include "Conventions/Constants.h"
#include "FluxDrivers/GAtmoFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
GAtmoFlux::GAtmoFlux()
{
  fInitialized = 0;
}
//___________________________________________________________________________
GAtmoFlux::~GAtmoFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
double GAtmoFlux::MaxEnergy(void)
{
  return TMath::Min(fMaxEv, fMaxEvCut);
}
//___________________________________________________________________________
bool GAtmoFlux::GenerateNext(void)
{
  while(1) {
     // Attempt to generate next flux neutrino
     bool nextok = this->GenerateNext_1try();
     if(!nextok) continue;

     // Check generated neutrino energy against max energy.
     // We may have to reject the current neutrino if a user-defined max
     // energy cut restricts the available range of energies.
     const TLorentzVector & p4 = this->Momentum();
     double E    = p4.Energy();
     double Emin = this->MinEnergy();
     double Emax = this->MaxEnergy();
     double wght = this->Weight();

     bool accept = (E<=Emax && E>=Emin && wght>0);
     if(accept) return true;
  }
  return false;
}
//___________________________________________________________________________
bool GAtmoFlux::GenerateNext_1try(void)
{
  // Must have run intitialization
  assert(fInitialized);

  // Reset previously generated neutrino code / 4-p / 4-x
  this->ResetSelection();

  // Get a RandomGen instance
  RandomGen * rnd = RandomGen::Instance();

  // Generate a (Ev, costheta) pair from the 'combined' flux histogram
  // and select (phi) uniformly over [0,2pi]
  double Ev       = 0.;
  double costheta = 0.;
  double phi      = 2.*kPi* rnd->RndFlux().Rndm();
  double weight   = 0;
  int    nu_pdg   = 0;

  if(fGenWeighted) {

     //
     // generate weighted flux
     //

     // generate events according to a power law spectrum,
     // then weight events by flux and inverse power law
     // (note: cannot use index alpha=1)
     double alpha = fSpectralIndex; 

     double emin = TMath::Power(fEnergyBins[0],1.0-alpha);
     double emax = TMath::Power(fEnergyBins[fNumEnergyBins],1.0-alpha);
     Ev          = TMath::Power(emin+(emax-emin)*rnd->RndFlux().Rndm(),1.0/(1.0-alpha));
     costheta    = -1+2*rnd->RndFlux().Rndm();

     unsigned int nnu = fPdgCList->size();
     unsigned int inu = rnd->RndFlux().Integer(nnu);
     nu_pdg   = (*fPdgCList)[inu];

     if(Ev < fEnergyBins[0]) {
        LOG("Flux", pINFO) << "E < Emin";
	return false;
     }
 
     double flux = GetFlux( nu_pdg, Ev, costheta );

     if(flux<=0) {
        LOG("Flux", pINFO) << "Flux <= 0";
	return false;
     }

     weight = flux*TMath::Power(Ev,alpha);
  } 
  else {

     //
     // generate nominal flux
     //

     Axis_t ax = 0, ay = 0;
     fFluxSum2D->GetRandom2(ax, ay);
     Ev       = (double)ax;
     costheta = (double)ay;
     nu_pdg   = this->SelectNeutrino(Ev, costheta);
     weight   = 1.0;
  }

  // Compute etc trigonometric numbers
  double sintheta  = TMath::Sqrt(1-costheta*costheta);
  double cosphi    = TMath::Cos(phi);
  double sinphi    = TMath::Sin(phi);

  // Set the neutrino pdg code
  fgPdgC = nu_pdg;

  // Set the neutrino weight
  fWeight = weight;

  // Compute the neutrino momentum
  // The `-1' means it is directed towards the detector.
  double pz = -1.* Ev * costheta;
  double py = -1.* Ev * sintheta * cosphi;
  double px = -1.* Ev * sintheta * sinphi;

  // Default vertex is at the origin
  double z = 0.0;
  double y = 0.0;
  double x = 0.0;

  // Shift the neutrino position onto the flux generation surface.
  // The position is computed at the surface of a sphere with R=fRl
  // at the topocentric horizontal (THZ) coordinate system.
  if( fRl>0.0 ){
    z += fRl * costheta;
    y += fRl * sintheta * cosphi;
    x += fRl * sintheta * sinphi;
  }

  // Apply user-defined rotation from THZ -> user-defined topocentric 
  // coordinate system.
  if( !fRotTHz2User.IsIdentity() )
  {
    TVector3 tx3(x, y, z );
    TVector3 tp3(px,py,pz);

    tx3 = fRotTHz2User * tx3;
    tp3 = fRotTHz2User * tp3;

    x  = tx3.X();
    y  = tx3.Y();
    z  = tx3.Z();
    px = tp3.X();
    py = tp3.Y();
    pz = tp3.Z();
  }

  // If the position is left as is, then all generated neutrinos
  // would point towards the origin.
  // Displace the position randomly on the surface that is
  // perpendicular to the selected point P(xo,yo,zo) on the sphere
  if( fRt>0.0 ){
    TVector3 vec(x,y,z);               // vector towards selected point
    TVector3 dvec1 = vec.Orthogonal(); // orthogonal vector
    TVector3 dvec2 = dvec1;            // second orthogonal vector
    dvec2.Rotate(-kPi/2.0,vec);        // rotate second vector by 90deg, 
                                       // now forming a new orthogonal cartesian coordinate system
    double psi = 2.*kPi* rnd->RndFlux().Rndm(); // rndm angle [0,2pi]
    double random = rnd->RndFlux().Rndm();      // rndm number  [0,1]
    dvec1.SetMag(TMath::Sqrt(random)*fRt*TMath::Cos(psi));
    dvec2.SetMag(TMath::Sqrt(random)*fRt*TMath::Sin(psi));
    x += dvec1.X() + dvec2.X();
    y += dvec1.Y() + dvec2.Y();
    z += dvec1.Z() + dvec2.Z();
  }

  // Set the neutrino momentum and position 4-vectors with values
  // calculated at previous steps.
  fgP4.SetPxPyPzE(px, py, pz, Ev);
  fgX4.SetXYZT   (x,  y,  z,  0.);

  // Increment flux neutrino counter used for sample normalization purposes.
  fNNeutrinos++;

  // Report and exit
  LOG("Flux", pINFO)
       << "Generated neutrino: "
       << "\n pdg-code: " << fgPdgC
       << "\n p4: " << utils::print::P4AsShortString(&fgP4)
       << "\n x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
long int GAtmoFlux::NFluxNeutrinos(void) const
{
  return fNNeutrinos;
}
//___________________________________________________________________________
void GAtmoFlux::ForceMinEnergy(double emin)
{
  emin = TMath::Max(0., emin);
  fMinEvCut = emin;
}
//___________________________________________________________________________
void GAtmoFlux::ForceMaxEnergy(double emax)
{
  emax = TMath::Max(0., emax);
  fMaxEvCut = emax;
}
//___________________________________________________________________________
void GAtmoFlux::Clear(Option_t * opt)
{
// Dummy clear method needed to conform to GFluxI interface 
//
  LOG("Flux", pERROR) << "No clear method implemented for option:"<< opt;
}
//___________________________________________________________________________
void GAtmoFlux::GenerateWeighted(bool gen_weighted)
{
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GAtmoFlux:: SetSpectralIndex(double index)
{
  if( index != 1.0 ){
    fSpectralIndex = index;
  }
  else {
    LOG("Flux", pWARN) << "Warning: cannot use a spectral index of unity";
  }

  LOG("Flux", pNOTICE) << "Using Spectral Index = " << index;
}
//___________________________________________________________________________
void GAtmoFlux::SetUserCoordSystem(TRotation & rotation)
{
  fRotTHz2User = rotation;
}
//___________________________________________________________________________
void GAtmoFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing atmospheric flux driver";

  bool allow_dup = false;
  fPdgCList = new PDGCodeList(allow_dup);

  // initializing flux TH2D histos [ flux = f(Ev,costheta) ] & files
  fFluxFile.clear();
  fFlux2D.clear();
  fFluxSum2D = 0;
  fFluxSum2DIntg = 0;

  // setting maximum energy in flux files
  assert(fEnergyBins);
  fMaxEv = fEnergyBins[fNumEnergyBins];

  // Default option is to generate unweighted flux neutrinos
  // (flux = f(E,costheta) will be used as PDFs)
  // User can enable option to generate weighted neutrinos
  // (neutrinos will be generated uniformly over costheta, 
  // and using a power law function in neutrino energy.
  // The input flux = f(E,costheta) will be used for calculating a weight).
  // Using a weighted flux avoids statistical fluctuations at high energies.
  fSpectralIndex = 2.0;

  // weighting switched off by default
  this->GenerateWeighted(false);

  // Default: No min/max energy cut
  this->ForceMinEnergy(0.);
  this->ForceMaxEnergy(9999999999.);

  // Default radii
  fRl = 0.0;
  fRt = 0.0;

  // Default detector coord system: Topocentric Horizontal Coordinate system
  fRotTHz2User.SetToIdentity(); 

  // Reset `current' selected flux neutrino
  this->ResetSelection();

  // Reset number of neutrinos thrown so far
  fNNeutrinos = 0;

  // Done!
  fInitialized = 1;
}
//___________________________________________________________________________
void GAtmoFlux::ResetSelection(void)
{
// initializing running neutrino pdg-code, 4-position, 4-momentum

  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);
}
//___________________________________________________________________________
void GAtmoFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  map<int,TH2D*>::iterator rawiter = fFluxRaw2D.begin();
  for( ; rawiter != fFluxRaw2D.end(); ++rawiter) {
    TH2D * flux_histogram = rawiter->second;
    if(flux_histogram) {
       delete flux_histogram;
       flux_histogram = 0;
    }
  }
  fFluxRaw2D.clear();

  map<int,TH2D*>::iterator iter = fFlux2D.begin();
  for( ; iter != fFlux2D.end(); ++iter) {
    TH2D * flux_histogram = iter->second;
    if(flux_histogram) {
       delete flux_histogram;
       flux_histogram = 0;
    }
  }
  fFlux2D.clear();

  if (fFluxSum2D) delete fFluxSum2D;
  if (fPdgCList ) delete fPdgCList;

  delete [] fCosThetaBins;
  delete [] fEnergyBins;
}
//___________________________________________________________________________
void GAtmoFlux::SetRadii(double Rlongitudinal, double Rtransverse)
{
  LOG ("Flux", pNOTICE) << "Setting R[longitudinal] = " << Rlongitudinal;
  LOG ("Flux", pNOTICE) << "Setting R[transverse]   = " << Rtransverse;

  fRl = Rlongitudinal;
  fRt = Rtransverse;
}
//___________________________________________________________________________
void GAtmoFlux::AddFluxFile(int nu_pdg, string filename)
{
  if ( pdg::IsNeutrino(nu_pdg) || pdg::IsAntiNeutrino(nu_pdg) ) {
    fFluxFlavour.push_back(nu_pdg); fFluxFile.push_back(filename);
  } else {
    LOG ("Flux", pWARN) 
     << "Input particle code: " << nu_pdg << " not a neutrino!";
  }
}
//___________________________________________________________________________
void GAtmoFlux::SetFluxFile(int nu_pdg, string filename)
{
  return AddFluxFile( nu_pdg, filename );
}
//___________________________________________________________________________
bool GAtmoFlux::LoadFluxData(void)
{
  LOG("Flux", pNOTICE)
        << "Loading atmospheric neutrino flux simulation data";

  fFlux2D.clear();
  fPdgCList->clear();

  bool loading_status = true;

  for( unsigned int n=0; n<fFluxFlavour.size(); n++ ){
    int nu_pdg      = fFluxFlavour.at(n);
    string filename = fFluxFile.at(n);
    string pname = PDGLibrary::Instance()->Find(nu_pdg)->GetName();

    LOG("Flux", pNOTICE)
        << "Loading data for: " << pname;
  
    TH2D* hist = 0;

    std::map<int,TH2D*>::iterator myMapEntry = fFluxRaw2D.find(nu_pdg);

    if( myMapEntry != fFluxRaw2D.end() ){
      hist = myMapEntry->second;
    }
    
    if( hist==0 ){
      hist = this->CreateFluxHisto2D(pname.c_str(), pname.c_str());
      fFluxRaw2D.insert( map<int,TH2D*>::value_type(nu_pdg,hist) );
    }

    bool loaded = this->FillFluxHisto2D(hist, filename);

    loading_status = loading_status && loaded;
  }

  if(loading_status) {

    map<int,TH2D*>::iterator hist_iter = fFluxRaw2D.begin();
    for ( ; hist_iter != fFluxRaw2D.end(); ++hist_iter) {
      int   nu_pdg = hist_iter->first;
      TH2D* hist   = hist_iter->second;

      TH2D* hnorm = this->CreateNormalisedFluxHisto2D( hist );
      fFlux2D.insert( map<int,TH2D*>::value_type(nu_pdg,hnorm) );
      fPdgCList->push_back(nu_pdg);
    }

    LOG("Flux", pNOTICE)
          << "Atmospheric neutrino flux simulation data loaded!";
    this->AddAllFluxes();
    return true;
  }

  LOG("Flux", pERROR)
    << "Error loading atmospheric neutrino flux simulation data";
  return false;
}
//___________________________________________________________________________
TH2D* GAtmoFlux::CreateNormalisedFluxHisto2D(TH2D* h2)
{
  // sanity check
  if( h2==0 ) return 0;
  
  // make new histogram name
  TString histname = h2->GetName();
  histname.Append("_IntegratedFlux");  

  // make new histogram
  TH2D* hIntegratedFlux2D = (TH2D*)(h2->Clone(histname.Data()));
  hIntegratedFlux2D->Reset();

  // integrate flux in each bin
  Double_t dN_dEdS = 0.0;
  Double_t dS = 0.0;
  Double_t dE = 0.0;
  Double_t dN = 0.0;

  for( Int_t nx=0; nx<h2->GetXaxis()->GetNbins(); nx++ ){ // x-axis: energy
    for( Int_t ny=0; ny<h2->GetYaxis()->GetNbins(); ny++ ){ // y-axis: angle
      dN_dEdS = h2->GetBinContent(nx+1,ny+1);

      dE = h2->GetXaxis()->GetBinUpEdge(nx+1)
          - h2->GetXaxis()->GetBinLowEdge(nx+1);

      dS = 2.0*TMath::Pi()
         * ( h2->GetYaxis()->GetBinUpEdge(ny+1)
            - h2->GetYaxis()->GetBinLowEdge(ny+1) );

      dN = dN_dEdS*dE*dS;

      hIntegratedFlux2D->SetBinContent(nx+1,ny+1,dN);
    }
  }

  // return integrated flux
  return hIntegratedFlux2D; 
}
//___________________________________________________________________________
void GAtmoFlux::ZeroFluxHisto2D(TH2D * histo)
{
  LOG("Flux", pNOTICE) << "Forcing flux histogram contents to 0";

  for(unsigned int ie = 0; ie < fNumEnergyBins; ie++) {
    for(unsigned int ic = 0; ic < fNumCosThetaBins; ic++) {
       double energy   = fEnergyBins  [ie];
       double costheta = fCosThetaBins[ic];
       histo->Fill(energy,costheta,0.);
    }
  }
}
//___________________________________________________________________________
void GAtmoFlux::AddAllFluxes(void)
{
  LOG("Flux", pNOTICE)
       << "Computing combined flux & flux normalization factor";

  if(fFluxSum2D) delete fFluxSum2D;

  fFluxSum2D = this->CreateFluxHisto2D("sum", "combined flux" );

  map<int,TH2D*>::iterator iter = fFlux2D.begin();
  for( ; iter != fFlux2D.end(); ++iter) {
    TH2D * flux_histogram = iter->second;
    fFluxSum2D->Add(flux_histogram);
  }

  fFluxSum2DIntg = fFluxSum2D->Integral();
}
//___________________________________________________________________________
TH2D * GAtmoFlux::CreateFluxHisto2D(string name, string title)
{
  LOG("Flux", pNOTICE) << "Instantiating histogram: [" << name << "]";
  TH2D * h2 = new TH2D(
           name.c_str(), title.c_str(),
           fNumEnergyBins, fEnergyBins, fNumCosThetaBins, fCosThetaBins);
  return h2;
}
//___________________________________________________________________________
int GAtmoFlux::SelectNeutrino(double Ev, double costheta)
{
// Select a neutrino species at the input Ev and costheta given their
// relatve flux at this bin.
// Returns a neutrino PDG code

  unsigned int n = fPdgCList->size();
  double * flux = new double[n];

  unsigned int i=0;
  map<int,TH2D*>::iterator iter = fFlux2D.begin();
  for( ; iter != fFlux2D.end(); ++iter) {
     TH2D * flux_histogram = iter->second;
     int ibin = flux_histogram->FindBin(Ev,costheta);
     flux[i]  = flux_histogram->GetBinContent(ibin);
     i++;
  }
  double flux_sum = 0;
  for(i=0; i<n; i++) {
     flux_sum  += flux[i];
     flux[i]    = flux_sum;
     LOG("Flux", pDEBUG) 
       << "Sum{Flux(0->" << i <<")} = " << flux[i];
  }

  RandomGen * rnd = RandomGen::Instance();
  double R = flux_sum * rnd->RndFlux().Rndm();

  LOG("Flux", pDEBUG) << "R = " << R;

  for(i=0; i<n; i++) {
     if( R < flux[i] ) {
	delete [] flux;
        int selected_pdg = (*fPdgCList)[i];
        LOG("Flux", pINFO) 
          << "Selected neutrino PDG code = " << selected_pdg;
	return selected_pdg;
     }
  }

  // shouldn't be here
  LOG("Flux", pERROR) << "Could not select a neutrino species!";
  assert(false);

  return -1;
}

//___________________________________________________________________________
TH2D* GAtmoFlux::GetFluxHistogram(int flavour)
{
  TH2D* histogram = 0;

  std::map<int,TH2D*>::iterator myMapEntry = fFluxRaw2D.find(flavour);

  if( myMapEntry != fFluxRaw2D.end() ){
    histogram = myMapEntry->second;
  }

  return histogram;
} 

//___________________________________________________________________________
double GAtmoFlux::GetFlux(int flavour)
{
  TH2D* hFlux2D = (TH2D*)(GetFluxHistogram(flavour));

  if( hFlux2D==0 ) return 0.0;

  Double_t Flux = 0.0;
  Double_t dN_dEdS = 0.0;
  Double_t dS = 0.0;
  Double_t dE = 0.0;
  
  for( Int_t nx=0; nx<hFlux2D->GetXaxis()->GetNbins(); nx++ ){ // x-axis: energy
    for( Int_t ny=0; ny<hFlux2D->GetYaxis()->GetNbins(); ny++ ){ // y-axis: angle
      dN_dEdS = hFlux2D->GetBinContent(nx+1,ny+1);

      dE = hFlux2D->GetXaxis()->GetBinUpEdge(nx+1)
          - hFlux2D->GetXaxis()->GetBinLowEdge(nx+1);

      dS = 2.0*TMath::Pi()
         * ( hFlux2D->GetYaxis()->GetBinUpEdge(ny+1)
           - hFlux2D->GetYaxis()->GetBinLowEdge(ny+1) );

      Flux += dN_dEdS*dE*dS;
    }
  }

  return Flux;
}

//___________________________________________________________________________
double GAtmoFlux::GetFlux(int flavour, double energy)
{
  TH2D* hFlux2D = (TH2D*)(GetFluxHistogram(flavour));

  if( hFlux2D==0 ) return 0.0;

  Int_t nE = hFlux2D->GetXaxis()->FindBin(energy); 

  Double_t Flux = 0.0;
  Double_t dN_dEdS = 0.0;
  Double_t dS = 0.0;
  
  for( Int_t ny=0; ny<hFlux2D->GetYaxis()->GetNbins(); ny++ ){ // y-axis: angle
    dN_dEdS = hFlux2D->GetBinContent(nE,ny+1);

    dS = 2.0*TMath::Pi()
       * ( hFlux2D->GetYaxis()->GetBinUpEdge(ny+1)
         - hFlux2D->GetYaxis()->GetBinLowEdge(ny+1) );

    Flux += dN_dEdS*dS;
  }

  return Flux;
}

//___________________________________________________________________________
double GAtmoFlux::GetFlux(int flavour, double energy, double angle)
{
  TH2D* hFlux2D = (TH2D*)(GetFluxHistogram(flavour));

  if( hFlux2D==0 ) return 0.0; 

  Int_t nE = hFlux2D->GetXaxis()->FindBin(energy);
  Int_t nA = hFlux2D->GetYaxis()->FindBin(angle);
  
  return hFlux2D->GetBinContent(nE,nA);
}

