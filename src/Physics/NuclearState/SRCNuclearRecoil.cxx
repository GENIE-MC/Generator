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
#include "Physics/NuclearState/SRCNuclearRecoil.h"

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

//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil() :
SecondNucleonEmissionI("genie::SRCNuclearRecoil")
{

}
//___________________________________________________________________________
SRCNuclearRecoil::SRCNuclearRecoil(string config) :
  SecondNucleonEmissionI("genie::SRCNuclearRecoil", config )
{

}

//___________________________________________________________________________

SRCNuclearRecoil::~SRCNuclearRecoil()
{

}

//___________________________________________________________________________

void SRCNuclearRecoil::ProcessEventRecord(GHepRecord * evrec) const
{

  const Interaction *  interaction = evrec       -> Summary();
  const InitialState & init_state  = interaction -> InitState();
  const Target       & tgt         = init_state.Tgt();

  // do nothing for non-nuclear targets
  if(! tgt.IsNucleus()) return;

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  // Set this to either a proton or neutron to eject a secondary particle
  int eject_nucleon_pdg = this->SRCRecoilPDG( *nucleon, tgt );

  // Ejection of secondary particle
  if (eject_nucleon_pdg != 0) { EmitSecondNucleon(evrec,eject_nucleon_pdg); }

}

//___________________________________________________________________________

int SRCNuclearRecoil::SRCRecoilPDG( const GHepParticle & nucleon, const Target & tgt) const {

      int eject_nucleon_pdg = 0;

      // const int nucleus_pdgc = tgt->Pdg();
      const int nucleon_pdgc = nucleon.Pdg();

      double pN2 = TMath::Power(nucleon.P4()->Rho(),2.); // (momentum of struck nucleon)^2

      double kF = fNuclModel -> LocalFermiMomentum( tgt, 
						    nucleon_pdgc, 
						    nucleon.X4()->Vect().Mag() );

      if (TMath::Sqrt(pN2) > kF) {
        double Pp = (nucleon_pdgc == kPdgProton) ? fPPPairPercentage : fPNPairPercentage;
        RandomGen * rnd = RandomGen::Instance();
        double prob = rnd->RndGen().Rndm();
        eject_nucleon_pdg = (prob > Pp) ? kPdgNeutron : kPdgProton;
      }
      
      return eject_nucleon_pdg;
}
//___________________________________________________________________________
void SRCNuclearRecoil::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SRCNuclearRecoil::LoadConfig(void)
{

  SecondNucleonEmissionI::LoadConfig() ;

  this->GetParamDef("PNPairPercentage",       fPNPairPercentage,    0.95);

  if (fPNPairPercentage < 0. || fPNPairPercentage > 1.) { 

	LOG("SRCNuclearRecoil", pFATAL)
	<< "PNPairPercentage either less than 0 or greater than 1: Exiting" ;

	exit(78); 
  }

  fPPPairPercentage = 1. - fPNPairPercentage;

}
//____________________________________________________________________________
