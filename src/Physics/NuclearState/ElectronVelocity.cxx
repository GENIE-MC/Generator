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

#include <iostream>
#include <random>

using namespace genie;

//___________________________________________________________________________
ElectronVelocity::ElectronVelocity(const string & name, const string & config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
ElectronVelocity::~ElectronVelocity()
{

}
//___________________________________________________________________________
void ElectronVelocity::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a electron target
  if(!evrec->Summary()->ProcInfo().IsElectronScattering()) return;

  // give electron initial momentum
  this->InitializeVelocity(*evrec->Summary());

}
//___________________________________________________________________________
void ElectronVelocity::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ElectronVelocity::Configure(string config)
{
  std::cout<<config<<std::endl;
  Algorithm::Configure(config);
  this->LoadConfig();
}
// void ElectronVelocity::LoadConfig(void)
// {

// }

