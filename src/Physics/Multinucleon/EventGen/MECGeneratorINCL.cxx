//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Liang Liu <liangliu \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodeList.h"

#include "Physics/Multinucleon/EventGen/MECGeneratorINCL.h"
#include "Physics/NuclearState/NucleusGenI.h"


using namespace genie;
using namespace genie::controls;

//___________________________________________________________________________
MECGeneratorINCL::MECGeneratorINCL() :
MECGenerator("genie::MECGeneratorINCL", "Default")
{

}
//___________________________________________________________________________
MECGeneratorINCL::MECGeneratorINCL(string config) :
MECGenerator("genie::MECGeneratorINCL", config)
{

}
//___________________________________________________________________________
MECGeneratorINCL::~MECGeneratorINCL()
{

}
//___________________________________________________________________________
void MECGeneratorINCL::GenerateFermiMomentum(GHepRecord * event) const
{
// Generate the initial state di-nucleon cluster 4-momentum using the
// INCL-based nucleus generator, which samples position and momentum
// simultaneously and applies the configured h-potential (local-energy)
// correction internally.
//
  GHepParticle * target_nucleus  = event->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = event->HitNucleon();
  assert(nucleon_cluster);
  GHepParticle * remnant_nucleus = event->RemnantNucleus();
  assert(remnant_nucleus);

  // The new interface of nucleus model:
  // it generates the momentum and position of a cluster simultaneously.
  fNucleusGen->GenerateCluster(event);

  // GenerateCluster has already set HitNucleon's P4; re-extract for clarity.
  TVector3 p3 = event->HitNucleon()->P4()->Vect();

  LOG("FermiMover", pINFO)
     << "di-nucleon cluster momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();

  // target (initial) nucleus and nucleon cluster mass
  double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass(); // initial nucleus mass
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster->Pdg())-> Mass(); // nucleon cluster mass

  // nucleon cluster energy (on-shell)
  double EN = TMath::Sqrt(p3.Mag2() + M2n*M2n);

  // set the remnant nucleus and nucleon cluster 4-momenta at GHEP record
  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  TLorentzVector p4remnant  (-1*p3.Px(), -1*p3.Py(), -1*p3.Pz(),  Mi-EN);

  nucleon_cluster->SetMomentum(p4nclust);
  remnant_nucleus->SetMomentum(p4remnant);

  event->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(p4nclust);
}
//___________________________________________________________________________
void MECGeneratorINCL::GenerateNSVInitialHadrons(GHepRecord * event) const
{
// Same accept/reject envelope as MECGenerator::GenerateNSVInitialHadrons,
// but the di-nucleon cluster 4-momentum and the binding 4-momentum come
// from the INCL nucleus generator instead of two independent NuclearModelI
// nucleon draws plus scalar removal energies.

  LOG("MEC",pDEBUG) << "Generate Initial Hadrons - Start";

  Interaction * interaction = event->Summary();
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4nu(*neutrino->P4());

  // final state primary lepton & its 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // 4-momentum transfer
  TLorentzVector Q4 = p4nu - p4l;

  // target two-nucleon cluster and nucleus
  GHepParticle * target_nucleus = event->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * initial_nucleon_cluster = event->HitNucleon();
  assert(initial_nucleon_cluster);
  GHepParticle * remnant_nucleus = event->RemnantNucleus();
  assert(remnant_nucleus);

  // Work out the final di-nucleon cluster pdg from the initial pair + probe.
  int final_nucleon_cluster_pdg = 0;
  if (neutrino->Pdg() == 11) {
    // e-scattering: same cluster
    final_nucleon_cluster_pdg = initial_nucleon_cluster->Pdg();
  }
  else if (neutrino->Pdg() > 0) {
    if (initial_nucleon_cluster->Pdg() == kPdgClusterNP) {
      final_nucleon_cluster_pdg = kPdgClusterPP;
    }
    else if (initial_nucleon_cluster->Pdg() == kPdgClusterNN) {
      final_nucleon_cluster_pdg = kPdgClusterNP;
    }
    else {
      LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction "
                         << initial_nucleon_cluster->Pdg();
    }
  }
  else if (neutrino->Pdg() < 0) {
    if (initial_nucleon_cluster->Pdg() == kPdgClusterNP) {
      final_nucleon_cluster_pdg = kPdgClusterNN;
    }
    else if (initial_nucleon_cluster->Pdg() == kPdgClusterPP) {
      final_nucleon_cluster_pdg = kPdgClusterNP;
    }
    else {
      LOG("MEC", pERROR) << "Wrong pdg for a CC anti-neutrino MEC interaction "
                         << initial_nucleon_cluster->Pdg();
    }
  }

  TLorentzVector p4initial_cluster;
  TLorentzVector p4final_cluster;
  TLorentzVector tLVebind;

  bool accept = false;
  unsigned int iter = 0;

  while (!accept) {
    iter++;
    if (iter > kRjMaxIterations*1000) {
      LOG("MEC", pWARN)
        << "Couldn't select a valid W, Q^2 pair after "
        << iter << " iterations";
      event->EventFlags()->SetBitNumber(kKineGenErr, true);
      genie::exceptions::EVGThreadException exception;
      exception.SetReason("Couldn't select initial hadron kinematics");
      exception.SwitchOnFastForward();
      throw exception;
    }

    // Generate a di-nucleon cluster: position + correlated momentum, with
    // the configured INCL h-potential applied.
    fNucleusGen->GenerateCluster(event);

    p4initial_cluster = *(event->HitNucleon()->P4());

    // Binding 4-momentum from the nucleus generator (replaces the legacy
    // scalar removal-energy sum).
    tLVebind = fNucleusGen->GetClusterBindP4();

    // tLVebind is (or should be) negative on the energy axis.
    p4final_cluster = p4initial_cluster + Q4 + tLVebind;

    // Accept only if we can put the final cluster on shell.
    if (p4final_cluster.M() <
        PDGLibrary::Instance()->Find(final_nucleon_cluster_pdg)->Mass()) {
      accept = false;
    } else {
      accept = true;
    }
  }  // end accept loop

  // Write everything to the GHEP record.

  // Initial-state nucleon cluster
  initial_nucleon_cluster->SetMomentum(p4initial_cluster);

  // Remnant nucleus
  double Mi = PDGLibrary::Instance()->Find(target_nucleus->Pdg())->Mass();
  remnant_nucleus->SetMomentum(-1.0*p4initial_cluster.Px(),
                               -1.0*p4initial_cluster.Py(),
                               -1.0*p4initial_cluster.Pz(),
                               Mi - p4initial_cluster.E() - tLVebind.E());

  // Final (undecayed) di-nucleon cluster, at the interaction vertex.
  TLorentzVector v4(*neutrino->X4());

  GHepParticle p1(final_nucleon_cluster_pdg, kIStDecayedState,
                  2, -1, -1, -1, p4final_cluster, v4);

  event->AddParticle(p1);

  interaction->KinePtr()->SetHadSystP4(p4final_cluster);
}
//___________________________________________________________________________
void MECGeneratorINCL::LoadConfig(void)
{
  // Pull in everything the base class needs (fQ3Max, fSafetyFactor, the
  // SuSAv2 tolerance, and -- importantly -- the NuclearModel sub-algorithm,
  // which MECGenerator::LoadConfig asserts on).
  //
  // Even though the INCL path never calls fNuclModel (both methods that
  // would use it are overridden here), we still satisfy the base assert by
  // requiring the NuclearModel entry in the corresponding XML. This keeps
  // MECGenerator unchanged.
  MECGenerator::LoadConfig();

  fNucleusGen = nullptr;
  RgKey nuclgenkey = "NucleusGen";
  fNucleusGen = dynamic_cast<const NucleusGenI *>(this->SubAlg(nuclgenkey));
  assert(fNucleusGen);
}
//___________________________________________________________________________

#endif // __GENIE_INCL_ENABLED__
