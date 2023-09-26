///____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 \brief  It visits the event record & computes a Bohr Velocity for
          initial state electrons bound in coloumb potential.

 \author   Brinden Carlson <bcarlson1 \at ufl.edu>
          University of Florida & Fermilab
  
  \created December 5, 2022

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/ElectronVelocity.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/BohrElectronVelocity.h"
#include "Physics/NuclearState/ElectronVelocity.h"

#include <iostream>
#include <random>

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
BohrElectronVelocity::~BohrElectronVelocity()
{

}
BohrElectronVelocity::BohrElectronVelocity(const string & config) :
ElectronVelocity::ElectronVelocity("genie::BohrElectronVelocity", config)
{

}
BohrElectronVelocity::BohrElectronVelocity() : 
ElectronVelocity::ElectronVelocity("genie::BohrElectronVelocity","Default")
{

}

//___________________________________________________________________________
void BohrElectronVelocity::InitializeVelocity(Interaction & interaction) const{
  InitialState * init_state  = interaction.InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  unsigned int fZ = tgt->Z(); //Get Z value
  double nucleus_mass = tgt->Mass(); 
  double mu = reduced_mass(nucleus_mass); //Reduced mass
  double v = random_bohr_velocity(fZ,mu); //bohr velocity
  TVector3 v3 = randomize_direction_sphere(v); //Get spherically uniform random direciton
  double gamma = 1/sqrt(1-v3.Mag2()); //Get boost
  //Set 3 momentum
  auto p3 = kElectronMass*gamma*v3;
  //-- update the electron 4p
  TLorentzVector * p4 = tgt->HitPartP4Ptr(); //Initialize 4 momentum pointer
  p4->SetVectM( p3, kElectronMass);
}
//___________________________________________________________________________

//____________________________________________________________________________
double BohrElectronVelocity::reduced_mass(double nucleus_mass) const{
  //Return reduced mass of electron
  return kElectronMass*nucleus_mass/(kElectronMass+nucleus_mass);
}

unsigned BohrElectronVelocity::random_n(unsigned int fZ) const{
  //Isotope minimum z is 118 in accordance with maximum possible in GENIE
  std::array<int, 6> fnprobs {2,10,28,60,110,118}; //Cumulative Probability dist.

  //Get random electron orbital from atomic number Z
  RandomGen * rnd = RandomGen::Instance(); //Load seed 

  for (unsigned int i = 0; i<fnprobs.size(); i++){
    if (fZ < fnprobs[i]){
      fnprobs[i] = fZ; //Set orbital value to Z
    }
  }

  double x = fZ * rnd->RndDec().Rndm();
  unsigned int n = 0;
  unsigned int sel_n = 0;
  do {
    sel_n = n;
  } while (x > fnprobs[n++]);

  return sel_n+1; //Plus 1 to get to n
}

double BohrElectronVelocity::bohr_velocity(unsigned int fn, unsigned int fZ, double mu) const
{
  return fZ*kAem/fn * kElectronMass/mu;
}

double BohrElectronVelocity::random_bohr_velocity(unsigned int fZ, double mu) const{
  //Get random bohr velocity from n distribution
  unsigned int fn = random_n(fZ);
  return bohr_velocity(fn,fZ,mu);
}

TVector3 BohrElectronVelocity::randomize_direction_sphere(double v) const{
  RandomGen * rnd = RandomGen::Instance(); //Load seed
  TRandom3 gen = rnd->RndGen(); //Get generator
  double costheta = gen.Uniform(-1,1); //Polar
  double phi = gen.Uniform(0,2*kPi); //Azimuthal
  double sintheta = sqrt( 1 - costheta*costheta); // which is always positive because the angle is in the 0 -> pi range! it cannot be negative
  return TVector3(v*sintheta*TMath::Cos(phi),
                  v*sintheta*TMath::Sin(phi),
                  v*costheta); //Return vector with magnitude v and random direction
}


//___________________________________________________________________________