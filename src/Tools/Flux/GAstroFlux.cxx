//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 27, 2010 - CA
   This class was first added in 2.7.1.
 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear and GFluxI::Index as 
   these methods needed for pre-generation of flux interaction probabilities 
   in GMCJDriver. 

*/
//____________________________________________________________________________

#include <cassert>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GAstroFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
GAstroFlux::GAstroFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GAstroFlux::~GAstroFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
double GAstroFlux::MaxEnergy(void)
{
  return TMath::Min(kAstroDefMaxEv, fMaxEvCut);
}
//___________________________________________________________________________
bool GAstroFlux::GenerateNext(void)
{
  // Reset previously generated neutrino code / 4-p / 4-x
  this->ResetSelection();

  if(!fEnergySpectrum) {
    return false;
  }
  if(!fSolidAngleAcceptance) {
    return false;
  }
  if(fRelNuPopulations.size() == 0) {
    return false;
  }

  //
  // Generate neutrino energy & starting position at the Geocentric
  // coordinate system
  //

  double log10Emin = TMath::Log10(TMath::Max(kAstroDefMinEv,fMinEvCut));
  double log10Emax = TMath::Log10(TMath::Min(kAstroDefMaxEv,fMaxEvCut));

  double wght_species = 1.;
  double wght_energy  = 1.;
  double wght_origin  = 1.;

  int    nupdg     = 0;
  double log10E    = -99999;
  double phi       = -999999;
  double costheta  = -999999;

  bool status = true;

  status = fNuGen->SelectNuPdg(
     fGenWeighted, fRelNuPopulations, nupdg, wght_species);
  if(!status) {
     return false;
  }

  status = fNuGen->SelectEnergy(
     fGenWeighted, *fEnergySpectrum, log10Emin, log10Emax, log10E, wght_energy);
  if(!status) {
     return false;
  }
  double Ev = TMath::Power(10.,log10E);

  status = fNuGen->SelectOrigin(
    fGenWeighted, *fSolidAngleAcceptance, phi, costheta, wght_origin);
  if(!status) {
     return false;
  }

  //
  // Propagate through the Earth: Get position, 4-momentum and neutrino
  // pdg code at the boundary of the detector volume
  //

  status = fNuPropg->Go(phi, costheta, fDetCenter, fDetSize, nupdg, Ev);
  if(!status) {
     return false;
  }

  int        pnupdg = fNuPropg->NuPdgAtDetVolBoundary();
  TVector3 & px3    = fNuPropg->X3AtDetVolBoundary(); 
  TVector3 & pp3    = fNuPropg->P3AtDetVolBoundary(); 

  //
  // Rotate vectors: 

  // GEF translated to detector centre -> THZ
  px3 = fRotGEF2THz * px3;
  pp3 = fRotGEF2THz * pp3;

  // THZ -> Topocentric user-defined detetor system
  px3 = fRotTHz2User * px3;
  pp3 = fRotTHz2User * pp3;

  //
  // Set position, momentum, pdg code and weight variables reported back
  //
  fgWeight = wght_species * wght_energy * wght_origin;
  fgPdgC   = pnupdg;
  fgX4.SetVect(px3*(units::m/units::km));
  fgX4.SetT(0.);
  fgP4.SetVect(pp3);
  fgP4.SetE(pp3.Mag());

  return true;
}
//___________________________________________________________________________
void GAstroFlux::ForceMinEnergy(double emin)
{
  emin = TMath::Max(0., emin/units::GeV);
  fMinEvCut = emin;
}
//___________________________________________________________________________
void GAstroFlux::ForceMaxEnergy(double emax)
{
  emax = TMath::Max(0., emax/units::GeV);
  fMaxEvCut = emax;
}
//___________________________________________________________________________
void GAstroFlux::Clear(Option_t * opt)
{
// Dummy clear method needed to conform to GFluxI interface 
//
  LOG("Flux", pERROR) << "No clear method implemented for option:"<< opt;
}
//___________________________________________________________________________
void GAstroFlux::GenerateWeighted(bool gen_weighted)
{
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GAstroFlux::SetDetectorPosition(
      double latitude, double longitude, double depth, double size)
{
  depth = TMath::Max(0., depth/units::km);
  size  = TMath::Max(0., size /units::km);

  assert(latitude  >= -kPi/2. && latitude  <= kPi/2.);
  assert(longitude >= 0.      && longitude < 2.*kPi );

  // set inputs
  fDetGeoLatitude   = latitude;
  fDetGeoLongitude  = longitude;
  fDetGeoDepth      = depth; 
  fDetSize          = size;

  //
  // Compute detector/topocentric coordinate system center in the 
  // geocentric coordinate system.
  //

  double REarth = constants::kREarth/units::km;
  double radius = REarth - fDetGeoDepth;

  double theta = kPi/2. - fDetGeoLatitude;
  double phi   = fDetGeoLongitude;

  double sintheta = TMath::Sin(theta);
  double costheta = TMath::Cos(theta);
  double sinphi   = TMath::Sin(phi);
  double cosphi   = TMath::Cos(phi);

  double xdc = radius * sintheta * cosphi;
  double ydc = radius * sintheta * sinphi;
  double zdc = radius * costheta;

  fDetCenter.SetXYZ(xdc,ydc,zdc);

  //
  // Coordinate System Rotation:
  // GEF translated to detector centre -> THZ
  //
  // ...
  // ... TODO
  // ...
}
//___________________________________________________________________________
void GAstroFlux::SetRelNuPopulations(
  double nnue,    double nnumu,    double nnutau, 
  double nnuebar, double nnumubar, double nnutaubar)
{
  fRelNuPopulations.clear();
  fPdgCList->clear();

  if(nnue>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgNuE, nnue));
    fPdgCList->push_back(kPdgNuE);
  }
  if(nnumu>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgNuMu, nnumu));
    fPdgCList->push_back(kPdgNuMu);
  }
  if(nnutau>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgNuTau, nnutau));
    fPdgCList->push_back(kPdgNuTau);
  }
  if(nnuebar>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgAntiNuE, nnuebar));
    fPdgCList->push_back(kPdgAntiNuE);
  }
  if(nnumubar>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgAntiNuMu, nnumubar));
    fPdgCList->push_back(kPdgAntiNuMu);
  }
  if(nnutaubar>0.) {
    fRelNuPopulations.insert(
       map<int,double>::value_type(kPdgAntiNuTau, nnutaubar));
    fPdgCList->push_back(kPdgAntiNuTau);
  }

  double tot = nnue + nnumu + nnutau + nnuebar + nnumubar + nnutaubar;
  assert(tot>0.);

  // normalize to 1.
  map<int,double>::iterator iter = fRelNuPopulations.begin();
  for ( ; iter != fRelNuPopulations.end(); ++iter) {
    double fraction = iter->second;
    double norm_fraction = fraction/tot;
    fRelNuPopulations[iter->first] = norm_fraction;
  }
}
//___________________________________________________________________________
void GAstroFlux::SetEnergyPowLawIdx(double n)
{
  if(fEnergySpectrum) delete fEnergySpectrum;

  double log10Emin = TMath::Log10(kAstroDefMinEv);
  double log10Emax = TMath::Log10(kAstroDefMaxEv);

  fEnergySpectrum = 
    new TH1D("fEnergySpectrum","",kAstroNlog10EvBins,log10Emin,log10Emax);
  fEnergySpectrum->SetDirectory(0);

  for(int i=1; i<=kAstroNlog10EvBins; i++) {
    double log10E = fEnergySpectrum->GetBinCenter(i);
    double Ev     = TMath::Power(10., log10E);
    double flux   = TMath::Power(Ev, -1*n);
    fEnergySpectrum->SetBinContent(i,flux);
  }

  // normalize
  double max = fEnergySpectrum->GetMaximum();
  fEnergySpectrum->Scale(1./max);
}
//___________________________________________________________________________
void GAstroFlux::SetUserCoordSystem(TRotation & rotation)
{
  fRotTHz2User = rotation;
}
//___________________________________________________________________________
void GAstroFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing flux driver";

  bool allow_dup = false;
  fPdgCList = new PDGCodeList(allow_dup);

  // Default: No min/max energy cut
  this->ForceMinEnergy(kAstroDefMinEv);
  this->ForceMaxEnergy(kAstroDefMaxEv);

  // Generate weighted or un-weighted flux?
  fGenWeighted = true;

  // Detector position & size
  fDetGeoLatitude  = -1.; 
  fDetGeoLongitude = -1.; 
  fDetGeoDepth     = -1.; 
  fDetSize         = -1.;  
  fDetCenter.SetXYZ(0,0,0); // in the geocentric coord system

  // Normalized 2-D histogram (phi,costheta): detector solid angle
  // as seen from different positions across the face of the Earth
  fSolidAngleAcceptance = 0;

  // Neutrino energy spectrum
  // To be set via SetEnergyPowLawIdx()
  // Can be trivially modified to accomodate different spectra
  fEnergySpectrum = 0;

  // Relative neutrino populations
  // To be set via SetRelNuPopulations()
  // Default nue:numu:nutau:nuebar:numubar:nutaubar = 1:2:0:1:2:0
  fRelNuPopulations.clear();

  // Rotations
  fRotGEF2THz.SetToIdentity();
  fRotTHz2User.SetToIdentity();

  // Utility objects for generating and propagating neutrinos
  fNuGen   = new NuGenerator();
  fNuPropg = new NuPropagator(1.0*units::km);

  // Reset `current' selected flux neutrino
  this->ResetSelection();
}
//___________________________________________________________________________
void GAstroFlux::ResetSelection(void)
{
// initializing running neutrino pdg-code, 4-position, 4-momentum

  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);
}
//___________________________________________________________________________
void GAstroFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  fRelNuPopulations.clear();
  fPdgCList->clear();

  delete fPdgCList;
  if(fEnergySpectrum) delete fEnergySpectrum;
  if(fSolidAngleAcceptance) delete fSolidAngleAcceptance;

  delete fNuGen;
  delete fNuPropg;
}
//___________________________________________________________________________

//
// GAstroFlux utility class method definitions
// ..........................................................................
//
//___________________________________________________________________________
bool GAstroFlux::NuGenerator::SelectNuPdg (
   bool weighted, const map<int,double> & nupdgpdf, 
   int & nupdg, double & wght)
{
// select neutrino species based on relative neutrino species populations
//
  nupdg = 0;
  wght  = 0;

  if(nupdgpdf.size() == 0) {
    return false;
  }

  RandomGen * rnd = RandomGen::Instance();

  // Generate weighted flux:
  //
  if(weighted) {
     unsigned int nnu = nupdgpdf.size();
     unsigned int inu = rnd->RndFlux().Integer(nnu);
     map<int,double>::const_iterator iter = nupdgpdf.begin();
     advance(iter,inu);
     nupdg = iter->first;
     wght  = iter->second;
  }
  // Generate un-weighted flux:
  //
  else {
     double xsum = 0.;
     double xrnd = rnd->RndFlux().Uniform();
     map<int,double>::const_iterator iter = nupdgpdf.begin();
     for( ; iter != nupdgpdf.end(); ++iter) {
       xsum += iter->second;
       if(xrnd < xsum) {
         nupdg = iter->first;
         break;
       }
     }
     wght = 1.;
  }

  if(nupdg==0) {
    return false;
  }

  return true;
}
//___________________________________________________________________________
bool GAstroFlux::NuGenerator::SelectEnergy(
  bool weighted, TH1D & log10Epdf, double log10Emin, double log10Emax, 
  double & log10E, double & wght)
{
// select neutrino energy
//

  log10E   = -9999999;
  wght     = 0;

  if(log10Emax <= log10Emin) {
    return false;
  }

  // Generate weighted flux:
  //
  if(weighted) {
     RandomGen * rnd = RandomGen::Instance();
     log10E  = log10Emin + (log10Emax-log10Emin) * rnd->RndFlux().Rndm();
     wght    = log10Epdf.GetBinContent(log10Epdf.FindBin(log10E));
  } 

  // Generate un-weighted flux:
  //
  else {
     do {
       log10E = log10Epdf.GetRandom();
     }
     while(log10E < log10Emin || log10E > log10Emax);
     wght = 1.;
  }

  return true;
}
//___________________________________________________________________________
bool GAstroFlux::NuGenerator::SelectOrigin(
  bool weighted, TH2D & opdf, 
  double & phi, double & costheta, double & wght)
{
  wght     = 0;
  costheta = -999999;
  phi      = -999999;

  // Generate weighted flux:
  //
  if(weighted) {
     RandomGen * rnd = RandomGen::Instance();
     phi      = 2.*kPi * rnd->RndFlux().Rndm();
     costheta = -1. + 2.*rnd->RndFlux().Rndm();
     wght     = opdf.GetBinContent(opdf.FindBin(phi,costheta));
  } 

  // Generate un-weighted flux:
  //
  else {
     opdf.GetRandom2(phi,costheta);
     wght = 1.;
  }

  return true;
}
//___________________________________________________________________________
bool GAstroFlux::NuPropagator::Go(
  double phi, double costheta, const TVector3 & detector_centre, 
  double detector_sz, int nu_pdg, double Ev)
{
  // initialize neutrino code
  fNuPdg = nu_pdg;

  //
  // initialize neutrino position vector 
  //
  double sintheta  = TMath::Sqrt(1-costheta*costheta);
  double cosphi    = TMath::Cos(phi);
  double sinphi    = TMath::Sin(phi);
  double REarth    = constants::kREarth/units::km;
  double zs = REarth * costheta;
  double ys = REarth * sintheta * cosphi;
  double xs = REarth * sintheta * sinphi;

  TVector3 start_position(xs,ys,zs);
  fX3 = start_position - detector_centre;

  //
  // initialize neutrino momentum 4-vector 
  //
  TVector3 direction_unit_vec = -1. * fX3.Unit();
  fP3 = Ev * direction_unit_vec;

  //
  // step through the Earth
  //

  LOG("Flux", pWARN) << "|dist|    = " << fX3.Mag();
  LOG("Flux", pWARN) << "|detsize| = " << detector_sz;

  while(1) {
    double currdist = fX3.Mag();
    if(currdist <= detector_sz - 0.1) break;

    double stepsz = (currdist-detector_sz > fStepSize) ? 
                         fStepSize : currdist-detector_sz;
    if(stepsz <= 0.) break;

//    LOG("Flux", pWARN) << "Stepping by... |dr| = " << stepsz;

    //
    // check earth density at current position, calculate interaction
    // probability over next step size, decide whether it interacts and
    // what happens if it does...
    //
    // ... todo ...

    fX3 += (stepsz * direction_unit_vec);
  }

  return true;
}
//___________________________________________________________________________

//
// GDiffuseAstroFlux concrete flux driver
// ..........................................................................
//
//___________________________________________________________________________
GDiffuseAstroFlux::GDiffuseAstroFlux() :
GAstroFlux()
{

}
//___________________________________________________________________________
GDiffuseAstroFlux::~GDiffuseAstroFlux() 
{

}
//___________________________________________________________________________

//
// GPointSourceAstroFlux concrete flux driver
// ..........................................................................
//
//___________________________________________________________________________
GPointSourceAstroFlux::GPointSourceAstroFlux() :
GAstroFlux()
{
  fPntSrcName.clear(); 
  fPntSrcRA.  clear(); 
  fPntSrcDec. clear();
  fPntSrcRelI.clear();

  fPntSrcTotI = 0;
}
//___________________________________________________________________________
GPointSourceAstroFlux::~GPointSourceAstroFlux() 
{

}
//___________________________________________________________________________
bool GPointSourceAstroFlux::GenerateNext(void)
{
  return true;
}
//___________________________________________________________________________
void GPointSourceAstroFlux::AddPointSource(
   string name, double ra, double dec, double rel_intensity)
{
  bool accept = 
         (ra  >= 0.      && ra  < 2.*kPi)  &&
         (dec >= -kPi/2. && dec <= kPi/2.) &&
         (rel_intensity > 0) &&
         (name.size() > 0);

  if(accept) {
     int id = fPntSrcName.size();

     fPntSrcName.insert( map<int, string>::value_type(id, name         ) );
     fPntSrcRA.  insert( map<int, double>::value_type(id, ra           ) );
     fPntSrcDec. insert( map<int, double>::value_type(id, dec          ) );
     fPntSrcRelI.insert( map<int, double>::value_type(id, rel_intensity) );

     fPntSrcTotI += rel_intensity;
  }
}  
//___________________________________________________________________________
bool GPointSourceAstroFlux::SelectSource(void)
{
  if(fPntSrcRelI.size() == 0) {
    return false;
  }

  if(fPntSrcTotI <= 0.) {
    return false;
  }

  unsigned int srcid = 0;
  double wght  = 0;

  // Generate weighted flux:
  //
  if(fGenWeighted) {
     RandomGen * rnd = RandomGen::Instance();
     unsigned int nsrc = fPntSrcName.size();
     srcid = rnd->RndFlux().Integer(nsrc);
     wght  = fPntSrcRelI[srcid] / fPntSrcTotI;
  }
  // Generate un-weighted flux:
  //
  else {
     RandomGen * rnd = RandomGen::Instance();
     double xsum = 0.;
     double xrnd = fPntSrcTotI * rnd->RndFlux().Uniform();
     map<int,double>::const_iterator piter = fPntSrcRelI.begin();
     for( ; piter != fPntSrcRelI.end(); ++piter) {
       xsum += piter->second;
       if(xrnd < xsum) {
         srcid = piter->first;
         break;
       }
     }
     wght = 1.;
  }

  //
  fSelSourceId  = srcid;

  //
  fgWeight *= wght;

  return true;
}
//___________________________________________________________________________

