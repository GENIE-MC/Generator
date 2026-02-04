//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC
*/
//____________________________________________________________________________

#include <cassert>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GPowerLawFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GPowerLawFlux,genie::flux::GPowerLawFlux)

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GPowerLawFlux::GPowerLawFlux() :
GFluxI()
{
  // default ctor for consistency with GFluxDriverFactory needs
  // up to user to call Initialize() to set energy and flavor(s)
}
//____________________________________________________________________________
GPowerLawFlux::GPowerLawFlux(double alpha, double emin, double emax, int pdg) :
GFluxI()
{
  this->Initialize(alpha,emin,emax,pdg);
}
//___________________________________________________________________________
GPowerLawFlux::GPowerLawFlux(
       double alpha, double emin, double emax, const map<int,double> & numap) :
GFluxI()
{
  this->Initialize(alpha,emin,emax,numap);
}
//___________________________________________________________________________
GPowerLawFlux::~GPowerLawFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GPowerLawFlux::GenerateNext(void)
{
  RandomGen * rnd = RandomGen::Instance();
  double p = fProbMax * rnd->RndFlux().Rndm();

  map<int,double>::const_iterator iter;
  for(iter = fProb.begin(); iter != fProb.end(); ++iter) {
     int    nupdgc = iter->first;
     double prob   = iter->second;
     if (p<prob) {
       fgPdgC = nupdgc;
       break;
     }
  }

  double Ev = 0.;
  if (fSpectralIndex==1) Ev = TMath::Exp(TMath::Log(fMinEv)+rnd->RndFlux().Rndm()*TMath::Log(fMaxEv/fMinEv));
  else {
    double pemin = TMath::Power(fMinEv,  1.-fSpectralIndex);
    double pemax = TMath::Power(fMaxEv,  1.-fSpectralIndex);
    Ev = TMath::Power(pemin+(pemax-pemin)*rnd->RndFlux().Rndm(),1./(1.-fSpectralIndex));
  }

  fgP4.SetPxPyPzE (0.,0.,Ev,Ev);

  LOG("Flux", pINFO)
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GPowerLawFlux::Clear(Option_t * opt)
{
// Dummy clear method needed to conform to GFluxI interface
//
  LOG("Flux", pERROR) <<
      "No Clear(Option_t * opt) method implemented for opt: "<< opt;
}
//___________________________________________________________________________
void GPowerLawFlux::GenerateWeighted(bool gen_weighted)
{
// Dummy implementation needed to conform to GFluxI interface
//
  LOG("Flux", pERROR) <<
      "No GenerateWeighted(bool gen_weighted) method implemented for " <<
      "gen_weighted: " << gen_weighted;
}
//___________________________________________________________________________
void GPowerLawFlux::Initialize(double alpha, double emin, double emax, int pdg)
{
  map<int,double> numap;
  numap.insert( map<int, double>::value_type(pdg, 1.) );

  this->Initialize(alpha,emin,emax,numap);
}
//___________________________________________________________________________
void GPowerLawFlux::Initialize(double alpha, double emin, double emax, const map<int,double> & numap)
{
  LOG("Flux", pNOTICE) << "Initializing GPowerLawFlux driver";

  fSpectralIndex = alpha;
  fMinEv         = emin;
  fMaxEv         = emax;

  LOG("Flux", pNOTICE) << "Spectral Index : " << fSpectralIndex;
  LOG("Flux", pNOTICE) << "Energy range : " << fMinEv << " , " << fMaxEv;

  fPdgCList = new PDGCodeList;
  fPdgCList->clear();

  fProbMax = 0;
  fProb.clear();

  map<int,double>::const_iterator iter;
  for(iter = numap.begin(); iter != numap.end(); ++iter) {
        int    nupdgc = iter->first;
        double nuwgt  = iter->second;

        fPdgCList->push_back(nupdgc);

        fProbMax+=nuwgt;
        fProb.insert(map<int, double>::value_type(nupdgc,fProbMax));
  }

  fgPdgC = 0;
  fgX4.SetXYZT (0.,0.,0.,0.);
}
//___________________________________________________________________________
void GPowerLawFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList) delete fPdgCList;
}
//___________________________________________________________________________
void GPowerLawFlux::SetDirectionCos(double dx, double dy, double dz)
{
  TVector3 dircos1 = TVector3(dx,dy,dz).Unit();
  LOG("Flux", pNOTICE) << "SetDirectionCos "
                       << utils::print::P3AsString(&dircos1);
  double E = fgP4.E();
  fgP4.SetVect(E*dircos1);

}
//___________________________________________________________________________
void GPowerLawFlux::SetRayOrigin(double x, double y, double z)
{
  TVector3 xyz(x,y,z);
  LOG("Flux", pNOTICE) << "SetRayOrigin "
                       << utils::print::Vec3AsString(&xyz);
  fgX4.SetVect(xyz);
}
//___________________________________________________________________________
void GPowerLawFlux::SetNuDirection(const TVector3 & direction)
{
  SetDirectionCos(direction.x(), direction.y(), direction.z());
}
//___________________________________________________________________________
void GPowerLawFlux::SetBeamSpot(const TVector3 & spot)
{
  SetRayOrigin(spot.x(), spot.y(), spot.z());
}
//___________________________________________________________________________
