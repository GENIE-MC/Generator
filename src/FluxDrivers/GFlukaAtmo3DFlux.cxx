//____________________________________________________________________________
/*!

\class   genie::GFlukaAtmo3DFlux

\brief   A flux driver for the FLUKA 3-D Atmospheric Neutrino Flux

\ref     Astrop.Phys.19 (2003) p.269; hep-ph/0207035; hep-ph/9907408
         Alfredo.Ferrari     <Alfredo.Ferrari@cern.ch>
         Paola.Sala          <Paola.Sala@cern.ch>
         Giuseppe Battistoni <Giuseppe.Battistoni@mi.infn.it>
         Teresa Montaruli    <Teresa.Montaruli@ba.infn.it>

         To be able to use this flux driver you will need to download the
         flux data from:  http://lxmi.mi.infn.it/~battist/neutrino.html

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 3, 2005 [during the most boring MINOS shift ever!]

*/
//____________________________________________________________________________

#include <fstream>

#include <TH2D.h>
#include <TMath.h>

#include "Conventions/Constants.h"
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGCodes.h"

using std::ifstream;
using std::ios;
using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
GFlukaAtmo3DFlux::GFlukaAtmo3DFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GFlukaAtmo3DFlux::~GFlukaAtmo3DFlux()
{
  for(unsigned int iflux = 0;
        iflux < kNNu; iflux++) { if (fFlux2D[iflux]) delete fFlux2D[iflux]; }
  if (fFluxSum2D) delete fFluxSum2D;
}
//___________________________________________________________________________
const PDGCodeList & GFlukaAtmo3DFlux::FluxParticles(void)
{
  return (*fPdgCList);
}
//___________________________________________________________________________
double GFlukaAtmo3DFlux::MaxEnergy(void)
{
  return fMaxEv;
}
//___________________________________________________________________________
bool GFlukaAtmo3DFlux::GenerateNext(void)
{
  //-- Reset previously generated neutrino code / 4-p / 4-x
  fgPdgC = 0;
  fgP4.SetPxPyPzE(0.,0.,0.,0.);
  fgX4.SetXYZT(0.,0.,0.,0.);

  //-- Get a RandomGen instance
  RandomGen * rnd = RandomGen::Instance();

  //-- Generate a (Ev, costheta, phi) triplet
  Axis_t aEv   = 0;
  Axis_t aCos8 = 0;
  fFluxSum2D->GetRandom2( aEv, aCos8);
  double Ev   = (double)aEv;
  double cos8 = (double)aCos8;
  double phi  = 2.*kPi*(rnd->Random2().Rndm());

  //-- Select a neutrino species
  int iflux = this->SelectNeutrino(Ev, cos8);
  assert(iflux>=0 && iflux<4);

  fgPdgC = (*fPdgCList)[iflux];

  //-- Compute the neutrino 4-p
  double sin8 = TMath::Sqrt(1-cos8*cos8);
  double pz = -1.* Ev * cos8;
  double py = -1.* Ev * sin8 * TMath::Cos(phi);
  double px = -1.* Ev * sin8 * TMath::Sin(phi);

  fgP4.SetPxPyPzE(px, py, pz, Ev);

  //-- Compute neutrino 4-x

  // compute its position at the surface of a sphere with R=fRl
  double z = fRl * cos8;
  double y = fRl * sin8 * TMath::Cos(phi);
  double x = fRl * sin8 * TMath::Sin(phi);

  // If the position is left as is, then all generated neutrinos
  // would point towards the origin.
  // Displace the position randomly on the surface that is
  // perpendicular to the selected point of the sphere.

  // Note: for the vector (xo,yo,zo), the perpendicular plane
  // a*z + b*y + c*x = 0 has a=zo,b=yo,c=xo.
  double az = z;
  double ay = y;
  double ax = x;
  // The displacement is governed by the input transverse scale
  double Rx = -fRt + fRt * rnd->Random2().Rndm();
  double Ry = -fRt + fRt * rnd->Random2().Rndm();
  // Compute new displaced position
  x = x + Rx;
  y = y + Ry;
  if(az!= 0) z = -(ay*y+ax*x)/az;
  else       z = 0;

  fgX4.SetXYZT(x,y,z,0.);

  return true;
}
//___________________________________________________________________________
int GFlukaAtmo3DFlux::PdgCode(void)
{
  LOG("Flux", pINFO)
         << "Loading GFluka 3-D Atmo. (Battistoni et al.) Simulation Data";
  return fgPdgC;
}
//___________________________________________________________________________
const TLorentzVector & GFlukaAtmo3DFlux::Momentum(void)
{
  LOG("Flux", pINFO)
         << "Loading GFluka 3-D Atmo. (Battistoni et al.) Simulation Data";
  return fgP4;
}
//___________________________________________________________________________
const TLorentzVector & GFlukaAtmo3DFlux::Position(void)
{
  LOG("Flux", pINFO)
         << "Loading GFluka 3-D Atmo. (Battistoni et al.) Simulation Data";
  return fgX4;
}
//___________________________________________________________________________
void GFlukaAtmo3DFlux::Initialize(void)
{
  // maximum energy in flux files & list of neutrinos
  fMaxEv    = kGFlk3DEv[kNGFlk3DEv-1];
  fPdgCList = new PDGCodeList(4);
  fPdgCList->push_back(kPdgNuMu);
  fPdgCList->push_back(kPdgNuMuBar);
  fPdgCList->push_back(kPdgNuE);
  fPdgCList->push_back(kPdgNuEBar);

  // Filenames for flux files [you need to get them from Battistoni's web
  // page - see class documentation]
  // You need to Set...() them before LoadFluxData()
  // You do not have to set all 4 if you need one flux component skipped
  // (but you do have to set at least one)
  fNSkipped  = 0;
  for (unsigned int iflux=0; iflux<kNNu; iflux++) fFluxFile[iflux] = "";

  // initializing flux TH2D histos [ flux = f(Ev,costheta) ]
  for (unsigned int iflux=0; iflux<kNNu; iflux++) fFlux2D[iflux] = 0;
  fFluxSum2D = 0;

  // initializing running neutrino pdg-code, 4-position, 4-momentum
  fgPdgC = 0;
  fgP4.SetPxPyPzE(0.,0.,0.,0.);
  fgX4.SetXYZT(0.,0.,0.,0.);
}
//___________________________________________________________________________
void GFlukaAtmo3DFlux::SetRadii(double Rlongitudinal, double Rtransverse)
{
  fRl = Rlongitudinal;
  fRt = Rtransverse;
}
//___________________________________________________________________________
bool GFlukaAtmo3DFlux::LoadFluxData(void)
{
  LOG("Flux", pINFO)
         << "Loading GFluka 3-D Atmo. (Battistoni et al.) Simulation Data";

  this->CreateFluxHisto2D(fFlux2D[0],"nu_mu",    "GFluka 3D flux: numu"   );
  this->CreateFluxHisto2D(fFlux2D[1],"nu_mu_bar","GFluka 3D flux: numubar");
  this->CreateFluxHisto2D(fFlux2D[2],"nu_e",     "GFluka 3D flux: nue"    );
  this->CreateFluxHisto2D(fFlux2D[3],"nu_e_bar", "GFluka 3D flux: nuebar" );

  bool loading_status = true;
  for(unsigned int iflux = 0; iflux < kNNu; iflux++) {
    bool loaded = this->FillFluxHisto2D(fFlux2D[iflux], fFluxFile[iflux]);
    loading_status = loading_status && loaded;
  }
  if(loading_status) {
     LOG("Flux", pINFO) << "GFluka 3-D Atmo. Simulation Data Loaded";
     this->AddAllFluxes();
     return true;
  }

  LOG("Flux", pERROR) << "Error in loading GFluka Simulation Data";
  return false;
}
//___________________________________________________________________________
bool GFlukaAtmo3DFlux::FillFluxHisto2D(TH2D * histo, string filename)
{
  LOG("Flux", pINFO) << "Reading flux data from: " << filename;

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
void GFlukaAtmo3DFlux::ZeroFluxHisto2D(TH2D * histo)
{
  LOG("Flux", pINFO) << "Setting flux histo. contents to 0 [skipping flux]";

  for(unsigned int ie = 0; ie < kNGFlk3DEv; ie++) {
    for(unsigned int ic = 0; ic < kNGFlk3DCos; ic++) {
       double energy   = kGFlk3DEv[ie]  - 1.E-4;
       double costheta = kGFlk3DCos[ic] - 1.E-4;
       histo->Fill(energy,costheta,0.);
    }
  }
}
//___________________________________________________________________________
void GFlukaAtmo3DFlux::AddAllFluxes(void)
{
  LOG("Flux", pINFO) << "Computing combined flux";

  this->CreateFluxHisto2D(fFluxSum2D, "sum", "combined flux");

  for(unsigned int iflux=0;
         iflux<kNNu; iflux++) fFluxSum2D->Add(fFlux2D[iflux]);
}
//___________________________________________________________________________
void GFlukaAtmo3DFlux::CreateFluxHisto2D(TH2D* h2, string name, string title)
{
  if (h2) delete h2;
  h2 = new TH2D(name.c_str(), title.c_str(),
                           kNGFlk3DEv-1,kGFlk3DEv, kNGFlk3DCos-1,kGFlk3DCos);
}
//___________________________________________________________________________
int GFlukaAtmo3DFlux::SelectNeutrino(double Ev, double costheta)
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
  double R = flux_sum * rnd->Random2().Rndm();

  LOG("Flux", pINFO) << "R = " << R;

  for(unsigned int iflux = 0; iflux < kNNu; iflux++) {
     if( R < flux[iflux] ) return iflux;
  }
  return -1;
}
//___________________________________________________________________________
