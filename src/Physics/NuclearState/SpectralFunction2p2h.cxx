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
#include "Physics/NuclearState/SpectralFunction2p2h.h"

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
SpectralFunction2p2h::SpectralFunction2p2h() :
SecondNucleonEmissionI("genie::SpectralFunction2p2h")
{

}
//___________________________________________________________________________
SpectralFunction2p2h::SpectralFunction2p2h(string config) :
  SecondNucleonEmissionI("genie::SpectralFunction2p2h", config)
{

}
//___________________________________________________________________________
SpectralFunction2p2h::~SpectralFunction2p2h()
{

}
//___________________________________________________________________________
void SpectralFunction2p2h::ProcessEventRecord(GHepRecord * evrec) const
{

    Interaction *  interaction = evrec       -> Summary();
    InitialState * init_state  = interaction -> InitStatePtr();
    Target *       tgt         = init_state  -> TgtPtr();

    if ( tgt -> A() <= 2 ) return ;
    if ( tgt -> Z() < 2 ) return ;

    FermiMoverInteractionType_t interaction_type = fNuclModel->GetFermiMoverInteractionType();
    
    if ( interaction_type == kFermiMoveEffectiveSF2p2h_eject ) {

        GHepParticle * nucleon = evrec->HitNucleon();
        int second_nucleon_pdg = nucleon->Pdg() == kPdgProton ? kPdgNeutron : kPdgProton ;
        SecondNucleonEmissionI::EmitSecondNucleon( evrec, second_nucleon_pdg );

    }

}
//____________________________________________________________________________
void SpectralFunction2p2h::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunction2p2h::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunction2p2h::LoadConfig(void)
{
  SecondNucleonEmissionI::LoadConfig() ;
}
//____________________________________________________________________________
