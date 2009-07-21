//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 04, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 08, 2008 - CA
   This trivial case was added in 2.3.1 so that single energy neutrinos can
   be easily used with the event generation driver that can handle a 
   target mix or detailed geometries.
 @ June 2, 2008 - CA
   Fix bug in Initialize() where weight was used as int
*/
//____________________________________________________________________________

#include <cassert>

#include "Conventions/Constants.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GMonoEnergeticFlux::GMonoEnergeticFlux(double Ev, int pdg) :
GFluxI()
{
  map<int,double> numap;
  numap.insert( map<int, double>::value_type(pdg, 1.) );

  this->Initialize(Ev,numap);
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
