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
BohrElectronVelocity::BohrElectronVelocity() : ElectronVelocity::ElectronVelocity()
{

}

//___________________________________________________________________________
void BohrElectronVelocity::InitializeVelocity(Interaction & interaction) const{
  InitialState * init_state  = interaction.InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  //Get random generator from genie
  RandomGen * rnd = RandomGen::Instance();
  TRandom3 gen = rnd->RndGen();

  unsigned int fZ = tgt->Z(); //Get Z value
  TVector3 v3;
  gen.Sphere(v3[0],v3[1],v3[2],random_bohr_velocity(fZ)); //Randomize direction with sphere - radius = bohr velocity
  double gamma = 1/sqrt(1-v3.Mag2()); //Get boost
  //Set 3 momentum
  auto p3 = kElectronMass*gamma*v3;
  //-- update the electron 4p
  TLorentzVector * p4 = tgt->HitEleP4Ptr(); //Initialize 4 momentum pointer
  p4->SetVectM( p3, kElectronMass2);
}
//___________________________________________________________________________

//____________________________________________________________________________

unsigned BohrElectronVelocity::random_n(unsigned int fZ) const{
  //Isotope minimum z is 118 in accordance with maximum possible in GENIE
  static int fMaxElectrons = 118; //Max total electrons
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

double BohrElectronVelocity::bohr_velocity(unsigned int fn, unsigned int fZ) const
{
  return fZ*kAem/fn;
}

double BohrElectronVelocity::random_bohr_velocity(unsigned int fZ) const{
  //Get random bohr velocity from n distribution
  unsigned int fn = random_n(fZ);
  return bohr_velocity(fn,fZ);
}


//___________________________________________________________________________