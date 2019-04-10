//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Joe Johnston (Univ of Pittsburgh) added code (Mar 18, 2016) to use
         either local or relativistic Fermi Gas for Pauli blocking.

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/NuclearState/PauliBlocker.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
PauliBlocker::PauliBlocker() :
EventRecordVisitorI("genie::PauliBlocker")
{

}
//___________________________________________________________________________
PauliBlocker::PauliBlocker(string config) :
EventRecordVisitorI("genie::PauliBlocker",  config)
{

}
//___________________________________________________________________________
PauliBlocker::~PauliBlocker()
{

}
//___________________________________________________________________________
void PauliBlocker::ProcessEventRecord(GHepRecord * evrec) const
{
  // Return if the neutrino was not scatterred off a nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("PauliBlock", pINFO)
	    << "No nuclear target found - The Pauli Blocker exits";
    return;
  }

  // Handle only QEL for now...
  // (can also be dark matter elastic)
  Interaction * interaction = evrec->Summary();
  const ProcessInfo & proc = interaction->ProcInfo();
  if(!proc.IsQuasiElastic() && !proc.IsDarkMatterElastic()) {
    LOG("PauliBlock", pINFO) << "Not a QEL event - The Pauli Blocker exits";
    return;
  }

  // Get the particle representing the initial hit nucleon
  GHepParticle * hit = evrec->HitNucleon();
  assert(hit);

  // Get the particle representing the recoiling final nucleon
  GHepParticle * recoil = evrec->Particle(hit->FirstDaughter());
  assert(recoil);
  int nuc_pdgc = recoil->Pdg();

  const Target& tgt = interaction->InitState().Tgt();
  double radius = hit->X4()->Vect().Mag();
  double kf = this->GetFermiMomentum(tgt, nuc_pdgc, radius);

  LOG("PauliBlock", pINFO) << "KF = " << kf;

  // get the recoil momentum
  double p = recoil->P4()->P(); // |p| for the recoil nucleon
  LOG("PauliBlock", pINFO) << "Recoil nucleon |P| = " << p;

   // check for pauli blocking
  bool is_blocked = (p < kf);

  // if it is blocked, set & thow an exception
  if ( is_blocked ) {
    LOG("PauliBlock", pNOTICE)
      << " *** The generated event is Pauli-blocked ("
      << "|p_{nucleon}| = " << p << " GeV < Fermi momentum = " << kf << " GeV) ***";

    evrec->EventFlags()->SetBitNumber(kPauliBlock, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Pauli-blocked event");

    // Include dark matter elastic
    if(proc.IsQuasiElastic() || proc.IsDarkMatterElastic()) {
      // nuclear suppression taken into account at the QEL cross
      // section - should attempt to regenerate the event as QEL
      exception.SwitchOnStepBack();
      exception.SetReturnStep(0);
    } else {
      // end this event generation thread and start again at the
      // interaction selection step
      // - this is irrelevant for the time being as we only handle QEL-
      exception.SwitchOnFastForward();
    }
    throw exception;
  }
}
//___________________________________________________________________________
void PauliBlocker::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadModelType();
}
//___________________________________________________________________________
void PauliBlocker::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadModelType();
}
//___________________________________________________________________________
void PauliBlocker::LoadModelType(void){
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Create a nuclear model object to check the model type
  RgKey nuclkey = "NuclearModel";
  RgAlg nuclalg = gc->GetAlg(nuclkey);
  AlgFactory * algf = AlgFactory::Instance();
  const NuclearModelI* nuclModel =
    dynamic_cast<const NuclearModelI*>(
			     algf->GetAlgorithm(nuclalg.name,nuclalg.config));
  // Check if the model is a local Fermi gas
  fLFG = (nuclModel && nuclModel->ModelType(Target()) == kNucmLocalFermiGas);

  if ( !fLFG ) {
    // get the Fermi momentum table for relativistic Fermi gas
	GetParam( "FermiMomentumTable", fKFTableName ) ;

    fKFTable = 0;

    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    fKFTable = kftp->GetTable(fKFTableName);
    assert(fKFTable);
  }
}
//___________________________________________________________________________
double PauliBlocker::GetFermiMomentum(const Target& tgt, int pdg_Nf,
  double radius) const
{
  // Pauli blocking should only be applied for nucleons
  assert( pdg::IsProton(pdg_Nf) || pdg::IsNeutron(pdg_Nf) );
  double kF = 0.;
  if ( fLFG ) {
    int A = tgt.A();
    bool is_p = pdg::IsProton( pdg_Nf );
    int numNuc = (is_p) ? tgt.Z() : tgt.N();
    double hbarc = kLightSpeed * kPlankConstant / units::fermi;
    kF = TMath::Power(3 * kPi2 * numNuc *
      genie::utils::nuclear::Density(radius, A), 1.0/3.0) * hbarc;
  }
  else {
    kF = fKFTable->FindClosestKF(tgt.Pdg(), pdg_Nf);
  }

  return kF;
}
