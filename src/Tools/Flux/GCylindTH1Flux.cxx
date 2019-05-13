//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - July 04, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index and 
   GFluxI::GenerateWeighted methods needed for pre-generation of flux
   interaction probabilities in GMCJDriver.

*/
//____________________________________________________________________________

#include <cassert>
#include <algorithm>

#include <TH1D.h>
#include <TF1.h>
#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GCylindTH1Flux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GCylindTH1Flux,genie::flux::GCylindTH1Flux)

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GCylindTH1Flux::GCylindTH1Flux()
{
  this->Initialize();
}
//___________________________________________________________________________
GCylindTH1Flux::~GCylindTH1Flux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GCylindTH1Flux::GenerateNext(void)
{
  //-- Reset previously generated neutrino code / 4-p / 4-x
  this->ResetSelection();

  //-- Generate an energy from the 'combined' spectrum histogram
  //   and compute the momentum vector
  double Ev = (double) fTotSpectrum->GetRandom();

  TVector3 p3(*fDirVec); // momentum along the neutrino direction
  p3.SetMag(Ev);         // with |p|=Ev

  fgP4.SetPxPyPzE(p3.Px(), p3.Py(), p3.Pz(), Ev);

  //-- Select a neutrino species from the flux fractions at the
  //   selected energy
  fgPdgC = (*fPdgCList)[this->SelectNeutrino(Ev)];

  //-- Compute neutrino 4-x

  if(fRt <= 0) {
    fgX4.SetXYZT(0.,0.,0.,0.);
  } 
  else {
    // Create a vector (vec) that points to a random position at a disk
    // of radius Rt passing through the origin, perpendicular to the
    // input direction.
    TVector3 vec0(*fDirVec);           // vector along the input direction
    TVector3 vec = vec0.Orthogonal();  // orthogonal vector

    double psi = this->GeneratePhi();  // rndm angle [0,2pi]
    double Rt  = this->GenerateRt();   // rndm R [0,Rtransverse]

    vec.Rotate(psi,vec0); // rotate around original vector
    vec.SetMag(Rt);       // set new norm

    // Set the neutrino position as beam_spot + vec
    double x = fBeamSpot->X() + vec.X();
    double y = fBeamSpot->Y() + vec.Y();
    double z = fBeamSpot->Z() + vec.Z();

    fgX4.SetXYZT(x,y,z,0.);
  }

  LOG("Flux", pINFO) << "Generated neutrino pdg-code: " << fgPdgC;
  LOG("Flux", pINFO)
        << "Generated neutrino p4: " << utils::print::P4AsShortString(&fgP4);
  LOG("Flux", pINFO)
             << "Generated neutrino x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GCylindTH1Flux::Clear(Option_t * opt)
{
// Dummy clear method needed to conform to GFluxI interface 
//
  LOG("Flux", pERROR) << 
      "No Clear(Option_t * opt) method implemented for opt: "<< opt;
}
//___________________________________________________________________________
void GCylindTH1Flux::GenerateWeighted(bool gen_weighted)
{
// Dummy implementation needed to conform to GFluxI interface
//
  LOG("Flux", pERROR) << 
      "No GenerateWeighted(bool gen_weighted) method implemented for " << 
      "gen_weighted: " << gen_weighted;
}
//___________________________________________________________________________
void GCylindTH1Flux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GCylindTH1Flux driver";

  fMaxEv       = 0;
  fPdgCList    = new PDGCodeList;
  fTotSpectrum = 0;
  fDirVec      = 0;
  fBeamSpot    = 0;
  fRt          =-1;
  fRtDep       = 0;

  this->ResetSelection();
  this->SetRtDependence("x");
  //eg, other example: this->SetRtDependence("pow(x,2)");
}
//___________________________________________________________________________
void GCylindTH1Flux::ResetSelection(void)
{
// initializing running neutrino pdg-code, 4-position, 4-momentum
  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);
}
//___________________________________________________________________________
void GCylindTH1Flux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fDirVec     ) delete fDirVec;
  if (fBeamSpot   ) delete fBeamSpot;
  if (fPdgCList   ) delete fPdgCList;
  if (fTotSpectrum) delete fTotSpectrum;
  if (fRtDep      ) delete fRtDep;

  unsigned int nspectra = fSpectrum.size();
  for(unsigned int i = 0; i < nspectra; i++) {
     TH1D * spectrum = fSpectrum[i];
     delete spectrum;
     spectrum = 0;
  }
}
//___________________________________________________________________________
void GCylindTH1Flux::SetNuDirection(const TVector3 & direction)
{
  if(fDirVec) delete fDirVec;
  fDirVec = new TVector3(direction);
}
//___________________________________________________________________________
void GCylindTH1Flux::SetBeamSpot(const TVector3 & spot)
{
  if(fBeamSpot) delete fBeamSpot;
  fBeamSpot = new TVector3(spot);
}
//___________________________________________________________________________
void GCylindTH1Flux::SetTransverseRadius(double Rt)
{
  LOG ("Flux", pNOTICE) << "Setting R[transverse] = " << Rt;
  fRt = Rt;

  if(fRtDep) fRtDep->SetRange(0,Rt);
}
//___________________________________________________________________________
void GCylindTH1Flux::AddEnergySpectrum(int nu_pdgc, TH1D * spectrum)
{
  LOG("Flux", pNOTICE) << "Adding flux spectrum for pdg = " << nu_pdgc;

  fPdgCList->push_back(nu_pdgc);

  bool accepted = (count(fPdgCList->begin(),fPdgCList->end(),nu_pdgc) == 1);
  if(!accepted) {
     LOG ("Flux", pWARN)
            << "The pdg-code isn't recognized and the spectrum was ignored";
  } else {
     fSpectrum.push_back(spectrum);

     int    nb  = spectrum->GetNbinsX();
     Axis_t max = spectrum->GetBinLowEdge(nb)+spectrum->GetBinWidth(nb);
     fMaxEv = TMath::Max(fMaxEv, (double)max);

     LOG("Flux", pNOTICE) 
          << "Updating maximum energy of flux particles to: " << fMaxEv;

     this->AddAllFluxes(); // update combined flux
  }
}
//___________________________________________________________________________
void GCylindTH1Flux::SetRtDependence(string rdep)
{
// Set the (functional form of) Rt dependence as string, eg "x*x+sin(x)"
// You do not need to set this method. The default behaviour is to generate
// flux neutrinos uniformly over the area of the cylinder's cross section.

  if(fRtDep) delete fRtDep;

  fRtDep = new TF1("rdep", rdep.c_str(), 0,fRt);
}
//___________________________________________________________________________
void GCylindTH1Flux::AddAllFluxes(void)
{
  LOG("Flux", pNOTICE) << "Computing combined flux";

  if(fTotSpectrum) delete fTotSpectrum;

  vector<TH1D *>::const_iterator spectrum_iter;

  unsigned int inu=0;
  for(spectrum_iter = fSpectrum.begin();
                       spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
     TH1D * spectrum = *spectrum_iter;

     if(inu==0) { fTotSpectrum = new TH1D(*spectrum); }
     else       { fTotSpectrum->Add(spectrum);        }
     inu++;
  }
}
//___________________________________________________________________________
int GCylindTH1Flux::SelectNeutrino(double Ev)
{
  const unsigned int n = fPdgCList->size();
  double fraction[n];

  vector<TH1D *>::const_iterator spectrum_iter;

  unsigned int inu=0;
  for(spectrum_iter = fSpectrum.begin();
                       spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
     TH1D * spectrum = *spectrum_iter;
     fraction[inu++] = spectrum->GetBinContent(spectrum->FindBin(Ev));
  }

  double sum = 0;
  for(inu = 0; inu < n; inu++) {
     sum += fraction[inu];
     fraction[inu] = sum;
     LOG("Flux", pDEBUG) << "SUM-FRACTION(0->" << inu <<") = " << sum;
  }

  RandomGen * rnd = RandomGen::Instance();
  double R = sum * rnd->RndFlux().Rndm();

  LOG("Flux", pDEBUG) << "R e [0,SUM] = " << R;

  for(inu = 0; inu < n; inu++) {if ( R < fraction[inu] ) return inu;}

  LOG("Flux", pERROR) << "Could not select a neutrino species";
  assert(false);

  return -1;
}
//___________________________________________________________________________
double GCylindTH1Flux::GeneratePhi(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double phi = 2.*kPi * rnd->RndFlux().Rndm(); // [0,2pi]
  return phi;
}
//___________________________________________________________________________
double GCylindTH1Flux::GenerateRt(void) const
{
  double Rt = fRtDep->GetRandom(); // rndm R [0,Rtransverse]
  return Rt;
}
//___________________________________________________________________________
