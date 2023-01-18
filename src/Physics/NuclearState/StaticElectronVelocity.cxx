///____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 
 \brief    It visits the event record & initializes a static velocity for
          initial state electron.

\author   Brinden Carlson <bcarlson1 \at ufl.edu>
          University of Florida & Fermilab

\created  December 5, 2022

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
StaticElectronVelocity::~StaticElectronVelocity()
{

}
StaticElectronVelocity::StaticElectronVelocity(const string & config) :
ElectronVelocity::ElectronVelocity("genie::StaticElectronVelocity", config)
{

}
StaticElectronVelocity::StaticElectronVelocity() : ElectronVelocity::ElectronVelocity()
{

}

//___________________________________________________________________________
void StaticElectronVelocity::InitializeVelocity(Interaction & interaction) const{
  InitialState * init_state  = interaction.InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  TLorentzVector * p4 = tgt->HitEleP4Ptr(); //Initialize 4 momentum pointer
  //These should be initialized like this by just in case
  TVector3 p3;
  p4->SetVectM(p3, kElectronMass2);
}
//___________________________________________________________________________

//____________________________________________________________________________
