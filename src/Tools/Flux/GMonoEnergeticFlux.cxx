//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - July 04, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 08, 2008 - CA
   This trivial case was added in 2.3.1 so that single energy neutrinos can
   be easily used with the event generation driver that can handle a 
   target mix or detailed geometries.
 @ Jun 02, 2008 - CA
   Fix bug in Initialize() where weight was used as int
 @ Jul 21, 2009 - RH
   Allow a bit more flexibility by giving the user the option to set the
   neutrino ray's direction cosines and origin.
 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index and 
   GFluxI::GenerateWeighted methods needed for pre-generation of flux
   interaction probabilities in GMCJDriver.

*/
//____________________________________________________________________________

#include <cassert>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GMonoEnergeticFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GMonoEnergeticFlux,genie::flux::GMonoEnergeticFlux)

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GMonoEnergeticFlux::GMonoEnergeticFlux() :
GFluxI()
{
  // default ctor for consistency with GFluxDriverFactory needs
  // up to user to call Initialize() to set energy and flavor(s)
}
//____________________________________________________________________________
GMonoEnergeticFlux::GMonoEnergeticFlux(double Ev, int pdg) :
GFluxI()
{
  this->Initialize(Ev,pdg);
}
//___________________________________________________________________________
GMonoEnergeticFlux::GMonoEnergeticFlux(
       double Ev, const map<int,double> & numap) :
GFluxI()
{
  this->Initialize(Ev,numap);
}
//___________________________________________________________________________
GMonoEnergeticFlux::~GMonoEnergeticFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GMonoEnergeticFlux::GenerateNext(void)
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

  LOG("Flux", pINFO) 
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GMonoEnergeticFlux::Clear(Option_t * opt)
{
// Dummy clear method needed to conform to GFluxI interface 
//
  LOG("Flux", pERROR) <<
      "No Clear(Option_t * opt) method implemented for opt: "<< opt;
}
//___________________________________________________________________________
void GMonoEnergeticFlux::GenerateWeighted(bool gen_weighted)
{
// Dummy implementation needed to conform to GFluxI interface
//
  LOG("Flux", pERROR) <<
      "No GenerateWeighted(bool gen_weighted) method implemented for " <<
      "gen_weighted: " << gen_weighted;
}
//___________________________________________________________________________
void GMonoEnergeticFlux::Initialize(double Ev, int pdg) 
{
  map<int,double> numap;
  numap.insert( map<int, double>::value_type(pdg, 1.) );

  this->Initialize(Ev,numap);
}
//___________________________________________________________________________
void GMonoEnergeticFlux::Initialize(double Ev, const map<int,double> & numap) 
{
  LOG("Flux", pNOTICE) << "Initializing GMonoEnergeticFlux driver";

  fMaxEv = Ev + 0.05;

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
  fgP4.SetPxPyPzE (0.,0.,Ev,Ev);
  fgX4.SetXYZT    (0.,0.,0.,0.);
}
//___________________________________________________________________________
void GMonoEnergeticFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList) delete fPdgCList;
}
//___________________________________________________________________________
void GMonoEnergeticFlux::SetDirectionCos(double dx, double dy, double dz)
{
  TVector3 dircos1 = TVector3(dx,dy,dz).Unit();
  LOG("Flux", pNOTICE) << "SetDirectionCos " 
                       << utils::print::P3AsString(&dircos1);
  double E = fgP4.E();
  fgP4.SetVect(E*dircos1);

}
//___________________________________________________________________________
void GMonoEnergeticFlux::SetRayOrigin(double x, double y, double z)
{
  TVector3 xyz(x,y,z);
  LOG("Flux", pNOTICE) << "SetRayOrigin "
                       << utils::print::Vec3AsString(&xyz);
  fgX4.SetVect(xyz);
}
//___________________________________________________________________________
void GMonoEnergeticFlux::SetNuDirection(const TVector3 & direction)
{
  SetDirectionCos(direction.x(), direction.y(), direction.z());
}
//___________________________________________________________________________
void GMonoEnergeticFlux::SetBeamSpot(const TVector3 & spot)
{
  SetRayOrigin(spot.x(), spot.y(), spot.z());
}
//___________________________________________________________________________
