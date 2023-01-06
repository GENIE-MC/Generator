///____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 Brinden Carlson
 University of Florida - August, 2022
 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory - October 08, 2004

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
//using namespace std;

//___________________________________________________________________________
// BohrElectronVelocity::BohrElectronVelocity()
// {
//
// }
BohrElectronVelocity::~BohrElectronVelocity()
{

}
BohrElectronVelocity::BohrElectronVelocity(const string & config) :
ElectronVelocity::ElectronVelocity("genie::BohrElectronVelocity", config)
{

}
BohrElectronVelocity::BohrElectronVelocity()
{

}

//___________________________________________________________________________
void BohrElectronVelocity::InitializeVelocity(Interaction & interaction) const{
  InitialState * init_state  = interaction.InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  //Get random generator from genie
  RandomGen * rnd = RandomGen::Instance();
  TRandom3 gen = rnd->RndGen();

  int Z = tgt->Z(); //Get Z value
  TLorentzVector * p4 = tgt->HitEleP4Ptr(); //Initialize 4 momentum pointer
  Double_t vx = 0;
  Double_t vy = 0;
  Double_t vz = 0;
  gen.Sphere(vx,vy,vz,random_bohr_velocity(Z)); //Randomize direction with sphere - radius = bohr velocity
  TVector3 v3(vx,vy,vz);
  float gamma = 1/sqrt(1-v3.Mag2()); //Get boost
  //Set 3 momentum
  TVector3 p3;
  p3.SetX(kElectronMass*gamma*v3.X());
  p3.SetY(kElectronMass*gamma*v3.Y());
  p3.SetZ(kElectronMass*gamma*v3.Z());
  //Calculate energy
  double PF2 = p3.Mag2(); //Magnitude of momentum 2
  double EN = sqrt(PF2+kElectronMass2);
  //-- update the electron 4p
  p4->SetPx( p3.Px() );
  p4->SetPy( p3.Py() );
  p4->SetPz( p3.Pz() );
  p4->SetE ( EN      );
  //std::cout<<"******"<<p3.Px()<<"******"<<p3.Py()<<"******"<<p3.Pz()<<"******"<<EN<<std::endl;
  //std::cout<<*init_state<<std::endl;
}
//___________________________________________________________________________

//____________________________________________________________________________

int BohrElectronVelocity::random_n(int Z) const{
  int fMaxElectrons = 118; //Max total electrons
  std::array<int, 6> fnprobs {2,10,28,60,110,118}; //Cumulative Probability dist.
  assert(Z<=fMaxElectrons); //Atomic number can't be greater than max electrons

  //Get random electron orbital from atomic number Z
  RandomGen * rnd = RandomGen::Instance(); //Load seed 

  for (unsigned int i = 0; i<fnprobs.size(); i++){
    if (Z < fnprobs[i]){
      fnprobs[i] = Z; //Set orbital value to Z
    }
  }

  double x = Z * rnd->RndDec().Rndm();
  int n = 0;
  int sel_n = 0;
  do {
    sel_n = n;
  } while (x > fnprobs[n++]);

  return sel_n+1; //Plus 1 to get to n
}

float BohrElectronVelocity::bohr_velocity(int n, int Z) const
{
  return Z*kAem/n;
}

float BohrElectronVelocity::random_bohr_velocity(int Z) const{
  //Get random bohr velocity from n distribution
  int n = random_n(Z);
  return bohr_velocity(n,Z);
}


//___________________________________________________________________________