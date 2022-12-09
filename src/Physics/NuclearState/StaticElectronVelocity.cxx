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
#include "Physics/NuclearState/StaticElectronVelocity.h"
#include "Physics/NuclearState/ElectronVelocity.h"

#include <iostream>
#include <random>

using namespace genie;
using namespace genie::constants;
//using namespace std;

//___________________________________________________________________________
// StaticElectronVelocity::StaticElectronVelocity()
// {
//
// }
StaticElectronVelocity::~StaticElectronVelocity()
{

}
StaticElectronVelocity::StaticElectronVelocity(const string & config) :
ElectronVelocity::ElectronVelocity("genie::StaticElectronVelocity", config)
{

}
StaticElectronVelocity::StaticElectronVelocity()
{

}

//___________________________________________________________________________
void StaticElectronVelocity::InitializeVelocity(Interaction & interaction) const{
  InitialState * init_state  = interaction.InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  //Get random generator from genie
  RandomGen * rnd = RandomGen::Instance();
  TRandom3 gen = rnd->RndGen();

  int Z = tgt->Z(); //Get Z value
  TLorentzVector * p4 = tgt->HitEleP4Ptr(); //Initialize 4 momentum pointer
  //These should be initialized like this by just in case
  p4->SetPx(0);
  p4->SetPy(0);
  p4->SetPz(0);
  p4->SetE ( kElectronMass);
}
//___________________________________________________________________________

//____________________________________________________________________________
