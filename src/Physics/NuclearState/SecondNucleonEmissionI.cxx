///____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Afroditi Papadopoulou <apapadop \at mit.edu>
         Massachusetts Institute of Technology - October 04, 2019

 @ October 4, 2019 - Afroditi Papadopoulou (AP)
   Created this new module that controls the addition of the recoil nucleon in the event record
   and extracts its kinematics
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
#include "Physics/NuclearState/SecondNucleonEmissionI.h"

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

using namespace genie;
using namespace genie::constants;

SecondNucleonEmissionI::SecondNucleonEmissionI(string name ) :
EventRecordVisitorI( name )
{

}
//___________________________________________________________________________
SecondNucleonEmissionI::SecondNucleonEmissionI(string name, string config) :
EventRecordVisitorI( name, config)
{

}
//___________________________________________________________________________
SecondNucleonEmissionI::~SecondNucleonEmissionI()
{

}

//___________________________________________________________________________
bool SecondNucleonEmissionI::EmitSecondNucleon( GHepRecord * evrec, const int eject_nucleon_pdg ) const {

  LOG("SecondNucleonEmissionI", pINFO) << "Adding a recoil nucleon with PDG " << eject_nucleon_pdg ;

  GHepParticle * nucleon = evrec->HitNucleon();

  GHepStatus_t status = kIStHadronInTheNucleus;
  int imom = evrec->TargetNucleusPosition();

  //-- Has opposite momentum from the struck nucleon
  double vx = nucleon->Vx();
  double vy = nucleon->Vy();
  double vz = nucleon->Vz();
  double px = -1.* nucleon->Px();
  double py = -1.* nucleon->Py();
  double pz = -1.* nucleon->Pz();
  double M  = PDGLibrary::Instance()->Find(eject_nucleon_pdg)->Mass();
  double E  = TMath::Sqrt(px*px+py*py+pz*pz+M*M);

  evrec->AddParticle( eject_nucleon_pdg, status, imom, -1, -1, -1, px, py, pz, E, vx, vy, vz, 0 );

  return true ;
}
//____________________________________________________________________________
void SecondNucleonEmissionI::LoadConfig(void)
{

  RgKey nuclkey = "NuclearModel";
  fNuclModel = 0;
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

}
//____________________________________________________________________________
