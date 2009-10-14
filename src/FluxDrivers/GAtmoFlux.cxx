//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Feb 05, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   This class was added in 2.3.1 by code factored out from the concrete
   GFlukaAtmo3DFlux driver & the newer, largely similar, GBartolAtmoFlux
   driver.
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
bool GAtmoFlux::GenerateNext(void)
{
  //-- Reset previously generated neutrino code / 4-p / 4-x
  this->ResetSelection();

  //-- Get a RandomGen instance
  RandomGen * rnd = RandomGen::Instance();

  //-- Generate a (Ev, costheta) pair from the 'combined' flux histogram
  //   and select (phi) uniformly over [0,2pi]
  Axis_t aEv   = 0;
  Axis_t aCos8 = 0;
  fFluxSum2D->GetRandom2( aEv, aCos8);
  double Ev   = (double)aEv;
  double cos8 = (double)aCos8;
  double phi  = 2.*kPi* rnd->RndFlux().Rndm();

  //-- etc trigonometric numbers
  double sin8   = TMath::Sqrt(1-cos8*cos8);
  double cosphi = TMath::Cos(phi);
  double sinphi = TMath::Sin(phi);

  //-- Select a neutrino species from the flux fractions at the selected
  //   (Ev,costtheta) pair

  fgPdgC = (*fPdgCList)[this->SelectNeutrino(Ev, cos8)];

  //-- Compute the neutrino 4-p
  //   (the - means it is directed towards the detector)
  double pz = -1.* Ev * cos8;
  double py = -1.* Ev * sin8 * cosphi;
  double px = -1.* Ev * sin8 * sinphi;

  fgP4.SetPxPyPzE(px, py, pz, Ev);

  //-- Compute neutrino 4-x

  // compute its position at the surface of a sphere with R=fRl
  double z = fRl * cos8;
  double y = fRl * sin8 * cosphi;
  double x = fRl * sin8 * sinphi;

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

  // displace the original vector & set the neutrino 4-position
  x += dvec.X();
  y += dvec.Y();
  z += dvec.Z();

  fgX4.SetXYZT(x,y,z,0.);

  LOG("Flux", pINFO)
       << "Generated neutrino: "
       << "\n pdg-code: " << fgPdgC
       << "\n p4: " << utils::print::P4AsShortString(&fgP4)
       << "\n x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GAtmoFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing base atmospheric flux driver";

  // setting maximum energy in flux files
  fMaxEv = fkEvBins[fkNEvBins];

  // setting neutrino-types in flux files
  PDGCodeList::size_type nnu = 4;
  fPdgCList = new PDGCodeList(nnu);
  (*fPdgCList)[0]= kPdgNuMu;
  (*fPdgCList)[1]= kPdgAntiNuMu;
  (*fPdgCList)[2]= kPdgNuE;
  (*fPdgCList)[3]= kPdgAntiNuE;

  // Filenames for flux files [you need to get them from the web
  // - see class documentation]
  // You need to Set...() them before LoadFluxData()
  // You do not have to set all 4 if you need one flux component skipped
  // (but you do have to set at least one)
  fNSkipped  = 0;
  for (unsigned int iflux=0; iflux<kNNu; iflux++) fFluxFile[iflux] = "";

  // initializing flux TH2D histos [ flux = f(Ev,costheta) ]
  //fFlux2D = new TH2D[kNNu];
  for (unsigned int iflux=0; iflux<kNNu; iflux++) fFlux2D[iflux] = 0;
  fFluxSum2D = 0;

  this->ResetSelection();
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

  for(unsigned int iflux = 0;
        iflux < kNNu; iflux++) { if (fFlux2D[iflux]) delete fFlux2D[iflux]; }
  if (fFluxSum2D) delete fFluxSum2D;
  if (fPdgCList ) delete fPdgCList;
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
void GAtmoFlux::SetCosThetaBins(unsigned int nbins, const double * bins)
{
  fkNCosBins = nbins;
  fkCosBins  = bins;
}
//___________________________________________________________________________
void GAtmoFlux::SetEnergyBins(unsigned int nbins, const double * bins)
{
  fkNEvBins = nbins;
  fkEvBins  = bins;
}
//___________________________________________________________________________
bool GAtmoFlux::LoadFluxData(void)
{
  LOG("Flux", pNOTICE) << "Creating Flux = f(Ev,cos8z) 2-D histograms";

  fFlux2D[0] = this->CreateFluxHisto2D("numu",   "atmo flux: numu"   );
  fFlux2D[1] = this->CreateFluxHisto2D("numubar","atmo flux: numubar");
  fFlux2D[2] = this->CreateFluxHisto2D("nue",    "atmo flux: nue"    );
  fFlux2D[3] = this->CreateFluxHisto2D("nuebar", "atmo flux: nuebar" );

  LOG("Flux", pNOTICE)
      << "Loading atmospheric neutrino flux simulation data";

  bool loading_status = true;
  for(unsigned int iflux = 0; iflux < kNNu; iflux++) {
    bool loaded = this->FillFluxHisto2D(fFlux2D[iflux], fFluxFile[iflux]);
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

  for(unsigned int ie = 0; ie < fkNEvBins; ie++) {
    for(unsigned int ic = 0; ic < fkNCosBins; ic++) {
       double energy   = fkEvBins[ie]  - 1.E-4;
       double costheta = fkCosBins[ic] - 1.E-4;
       histo->Fill(energy,costheta,0.);
    }
  }
}
//___________________________________________________________________________
void GAtmoFlux::AddAllFluxes(void)
{
  LOG("Flux", pNOTICE)
       << "Computing combined flux & flux normalization factor";

  fFluxSum2D = this->CreateFluxHisto2D("sum", "combined flux" );
  for(unsigned int iflux=0; iflux<kNNu; iflux++) {
	fFluxSum2D->Add(fFlux2D[iflux]);
  }
  fWeight = fFluxSum2D->Integral(0, fkNEvBins, 0, fkNCosBins, "width");
}
//___________________________________________________________________________
TH2D * GAtmoFlux::CreateFluxHisto2D(string name, string title)
{
  LOG("Flux", pNOTICE) << "Instantiating histogram: [" << name << "]";
  TH2D * h2 = new TH2D(
           name.c_str(), title.c_str(),
                     fkNEvBins, fkEvBins, fkNCosBins, fkCosBins);
  return h2;
}
//___________________________________________________________________________
int GAtmoFlux::SelectNeutrino(double Ev, double costheta)
{
  double flux[kNNu];
  for(unsigned int iflux=0; iflux<kNNu; iflux++) {
     int ibin = fFlux2D[iflux]->FindBin(Ev,costheta);
     flux[iflux] = fFlux2D[iflux]->GetBinContent(ibin);
  }
  double flux_sum = 0;
  for(unsigned int iflux=0; iflux<kNNu; iflux++) {
     flux_sum   += flux[iflux];
     flux[iflux] = flux_sum;
     LOG("Flux", pINFO) << "SUM-FLUX(0->" << iflux <<") = " << flux[iflux];
  }

  RandomGen * rnd = RandomGen::Instance();
  double R = flux_sum * rnd->RndFlux().Rndm();

  LOG("Flux", pINFO) << "R = " << R;

  for(unsigned int iflux = 0; iflux < kNNu; iflux++) {
     if( R < flux[iflux] ) return iflux;
  }

  LOG("Flux", pERROR) << "Could not select a neutrino species";
  assert(false);

  return -1;
}
//___________________________________________________________________________


