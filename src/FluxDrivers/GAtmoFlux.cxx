//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

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

     double log10emin = TMath::Log10(fEnergyBins[0]);
     double log10emax = TMath::Log10(fEnergyBins[fNumEnergyBins]);
     double dlog10e   = log10emax - log10emin;
     Ev       = TMath::Power(10., log10emin + dlog10e * rnd->RndFlux().Rndm());

     costheta = -1+2*rnd->RndFlux().Rndm();

     unsigned int nnu = fPdgCList->size();
     unsigned int inu = rnd->RndFlux().Integer(nnu);
     nu_pdg   = (*fPdgCList)[inu];

     if(Ev < fEnergyBins[0]) {
        LOG("Flux", pINFO) << "E < Emin";
	return false;
     }

     // calculate weight
     map<int,TH2D*>::iterator iter = fFlux2D.find(nu_pdg);
     if(iter == fFlux2D.end()) {
        LOG("Flux", pERROR) << "Can't find flux histogram for selected neutrino";
	return false;
     }
     TH2D* flux_histo = iter->second;
     if(!flux_histo) {
        LOG("Flux", pERROR) << "Null flux histogram!";
	return false;
     }
     int E_bin        = flux_histo->GetXaxis()->FindBin(Ev);
     int costheta_bin = flux_histo->GetYaxis()->FindBin(costheta);
     double flux      = flux_histo->GetBinContent(E_bin, costheta_bin);
   //double dE        = flux_histo->GetXaxis()->GetBinWidth(E_bin);
   //double dcostheta = flux_histo->GetYaxis()->GetBinWidth(costheta_bin);
   //if(fFluxSum2DIntg <= 0) {
   //   LOG("Flux", pERROR) << "Null flux integral!";
   //   return false;
   //}
   //weight = flux*dE*dcostheta / fFluxSum2DIntg;
     weight = flux;
  } 
  else {

     //
     // generate un-weighted flux
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

  // Compute the neutrino position (on the flux generation surface)
  // The position is computed at the surface of a sphere with R=fRl
  // at the topocentric horizontal (THZ) coordinate system.
  double z = fRl * costheta;
  double y = fRl * sintheta * cosphi;
  double x = fRl * sintheta * sinphi;

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
  TVector3 vec(x,y,z);              // vector towards selected point
  TVector3 dvec = vec.Orthogonal(); // orthogonal vector
  double psi = 2.*kPi* rnd->RndFlux().Rndm(); // rndm angle [0,2pi]
  double Rt  = fRt* rnd->RndFlux().Rndm();    // rndm norm  [0,Rtransverse]
  dvec.Rotate(psi,vec); // rotate around original vector
  dvec.SetMag(Rt);      // set new norm
  // displace the original vector 
  x += dvec.X();
  y += dvec.Y();
  z += dvec.Z();

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
void GAtmoFlux::GenerateWeighted(bool gen_weighted)
{
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GAtmoFlux::SetUserCoordSystem(TRotation & rotation)
{
  fRotTHz2User = rotation;
}
//___________________________________________________________________________
void GAtmoFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) 
    << "Initializing atmospheric flux driver";

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
  // (neutrinos will be generated uniformly over costheta, logE and
  // the input flux = f(E,costheta) will be used for calculating a weight).
  // Using a weighted flux avoids statistical fluctuations at high energies.
  this->GenerateWeighted(false);

  // Default: No min/max energy cut
  this->ForceMinEnergy(0.);
  this->ForceMaxEnergy(9999999999.);

  // Default detector coord system: Topocentric Horizontal Coordinate system
  fRotTHz2User.SetToIdentity(); 

  // Reset `current' selected flux neutrino
  this->ResetSelection();

  // Reset number of neutrinos thrown so far
  fNNeutrinos = 0;
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
void GAtmoFlux::SetFluxFile(int nu_pdg, string filename)
{
  if ( pdg::IsNeutrino(nu_pdg) || pdg::IsAntiNeutrino(nu_pdg) ) {
    fFluxFile.insert(map<int,string>::value_type(nu_pdg,filename));
  } else {
    LOG ("Flux", pWARN) 
     << "Input particle code: " << nu_pdg << " not a neutrino!";
  }
}
//___________________________________________________________________________
bool GAtmoFlux::LoadFluxData(void)
{
  LOG("Flux", pNOTICE)
        << "Loading atmospheric neutrino flux simulation data";

  fFlux2D.clear();
  fPdgCList->clear();

  bool loading_status = true;
  map<int,string>::iterator file_iter = fFluxFile.begin();
  for ( ; file_iter != fFluxFile.end(); ++file_iter) {
    int    nu_pdg    = file_iter->first;
    string filename  = file_iter->second;
    string pname = PDGLibrary::Instance()->Find(nu_pdg)->GetName();

    LOG("Flux", pNOTICE)
        << "Loading data for: " << pname;

    TH2D * hst = this->CreateFluxHisto2D(pname.c_str(), pname.c_str());

    bool loaded = this->FillFluxHisto2D(hst, filename);
    if(loaded) {
       fPdgCList->push_back(nu_pdg);
       fFlux2D.insert(map<int,TH2D*>::value_type(nu_pdg, hst));
    }
    loading_status = loading_status && loaded;
  }
  if(loading_status) {
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

  fFluxSum2DIntg = fFluxSum2D->Integral("width");
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


