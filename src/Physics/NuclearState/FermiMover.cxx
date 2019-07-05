///____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 08, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2009 - CA
   Added option to simulate the recoil nucleon due to short range corelation.
   The recoil nucleon is added in case that the selected hit nucleon has a
   momentum selected from the NN corelation tail. The 2 nucleons are put
   back-to-back. For the time-being using hard-coded relative fractions for
   the nn, pp, np pairs.
   The code for adding the recoil nuclear target at the GHEP record was moved
   into this processing step.

 @ Mar 18, 2016- Joe Johnston (SD)
   Call GenerateNucleon() with a target and a radius, so the local Fermi
   gas model can access the radius.
   Added a check to see if a local Fermi gas model is being used. If so,
   use a local Fermi gas model when deciding whether to eject a recoil nucleon.
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
#include "Physics/NuclearState/FermiMover.h"

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
FermiMover::FermiMover() :
EventRecordVisitorI("genie::FermiMover")
{

}
//___________________________________________________________________________
FermiMover::FermiMover(string config) :
EventRecordVisitorI("genie::FermiMover", config)
{

}
//___________________________________________________________________________
FermiMover::~FermiMover()
{

}
//___________________________________________________________________________
void FermiMover::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;

  // skip if no hit nucleon is set
  if(! evrec->HitNucleon()) return;

  // give hit nucleon a Fermi momentum
  this->KickHitNucleon(evrec);

  // add a recoiled nucleus remnant
  this->AddTargetNucleusRemnant(evrec);
}
//___________________________________________________________________________
void FermiMover::KickHitNucleon(GHepRecord * evrec) const
{
  Interaction *  interaction = evrec       -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  TLorentzVector * p4 = tgt->HitNucP4Ptr();

  // do nothing if the struct nucleon 4-momentum was set (eg as part of the
  // initial state selection)
  if(p4->Px()>0 || p4->Py()>0 || p4->Pz()>0) return;

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  // generate a Fermi momentum & removal energy
  // call GenerateNucleon with a radius in case the model is LocalFGM
  double rad = nucleon->X4()->Vect().Mag();
  fNuclModel->GenerateNucleon(*tgt,rad);

  TVector3 p3 = fNuclModel->Momentum3();
  double w    = fNuclModel->RemovalEnergy();

  LOG("FermiMover", pINFO)
     << "Generated nucleon momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();
  LOG("FermiMover", pINFO)
     << "Generated nucleon removal energy: w = " << w;

  double pF2 = p3.Mag2(); // (fermi momentum)^2

  nucleon->SetRemovalEnergy(w);

  // struck nucleon energy:
  // two possible prescriptions depending on whether you want to force
  // the sruck nucleon to be on the mass-shell or not...

  double EN=0;
  // Set this to either a proton or neutron to eject a secondary particle
  int eject_nucleon_pdg = 0;
  FermiMoverInteractionType_t interaction_type =
      fNuclModel->GetFermiMoverInteractionType();
  if (interaction_type == kFermiMoveEffectiveSF1p1h) {
    EN = nucleon->Mass() - w -
           pF2 / (2 * (nucleus->Mass() - nucleon->Mass()));
  } else if (interaction_type == kFermiMoveEffectiveSF2p2h_eject ||
             interaction_type == kFermiMoveEffectiveSF2p2h_noeject) {
    TParticlePDG * other_nucleon;
    if(nucleon->Pdg() == kPdgProton) {
      other_nucleon = PDGLibrary::Instance()->Find(kPdgNeutron);
    }
    else {
      other_nucleon = PDGLibrary::Instance()->Find(kPdgProton);
    }
    TParticlePDG * deuteron = PDGLibrary::Instance()->Find(1000010020);
    EN = deuteron->Mass() - 2 * w -
          TMath::Sqrt(pF2 + other_nucleon->Mass() * other_nucleon->Mass());
    if (interaction_type == kFermiMoveEffectiveSF2p2h_eject &&
        deuteron->PdgCode() != nucleus->Pdg()) {
      // Eject a remnant nucleon for the 2p2h process. Don't eject other
      // nucleon if target is deuterium- need to leave behind a
      // valid remnant nucleus.
      eject_nucleon_pdg = other_nucleon->PdgCode();
    }
  // Do default Fermi Moving
  } else  {
    if (!fKeepNuclOnMassShell) {
      //-- compute A,Z for final state nucleus & get its PDG code
      int nucleon_pdgc = nucleon->Pdg();
      bool is_p  = pdg::IsProton(nucleon_pdgc);
      int Z = (is_p) ? nucleus->Z()-1 : nucleus->Z();
      int A = nucleus->A() - 1;

      TParticlePDG * fnucleus = 0;
      int ipdgc = pdg::IonPdgCode(A, Z);
      fnucleus = PDGLibrary::Instance()->Find(ipdgc);
      if(!fnucleus) {
        LOG("FermiMover", pFATAL)
              << "No particle with [A = " << A << ", Z = " << Z
              << ", pdgc = " << ipdgc << "] in PDGLibrary!";
        exit(1);
      }
      //-- compute the energy of the struck (off the mass-shell) nucleus

      double Mf  = fnucleus -> Mass(); // remnant nucleus mass
      double Mi  = nucleus  -> Mass(); // initial nucleus mass

      EN = Mi - TMath::Sqrt(pF2 + Mf*Mf);
    } else {
      double MN  = nucleon->Mass();
      double MN2 = TMath::Power(MN,2);
      EN = TMath::Sqrt(MN2+pF2);
    }
    if (fSRCRecoilNucleon) {
      const int nucleus_pdgc = nucleus->Pdg();
      const int nucleon_pdgc = nucleon->Pdg();

      // Calculate the Fermi momentum, using a local Fermi gas if the
      // nuclear model is LocalFGM, and RFG otherwise
      double kF;
      if(fNuclModel->ModelType(*tgt) == kNucmLocalFermiGas){
	assert(pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc));
	int A = tgt->A();
	bool is_p = pdg::IsProton(nucleon_pdgc);
	double numNuc = (is_p) ? (double)tgt->Z():(double)tgt->N();
	double radius = nucleon->X4()->Vect().Mag();
	double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
	kF= TMath::Power(3*kPi2*numNuc*
		  genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc;
      }else{
	FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
	const FermiMomentumTable * kft  = kftp->GetTable("Default");
	kF = kft->FindClosestKF(nucleus_pdgc, nucleon_pdgc);
      }
      if (TMath::Sqrt(pF2) > kF) {
        double Pp = (nucleon->Pdg() == kPdgProton) ? 0.05 : 0.95;
        RandomGen * rnd = RandomGen::Instance();
        double prob = rnd->RndGen().Rndm();
        eject_nucleon_pdg = (prob > Pp) ? kPdgNeutron : kPdgProton;
      }
    }
  }

  //-- update the struck nucleon 4p at the interaction summary and at
  //   the GHEP record
  p4->SetPx( p3.Px() );
  p4->SetPy( p3.Py() );
  p4->SetPz( p3.Pz() );
  p4->SetE ( EN      );

  nucleon->SetMomentum(*p4); // update GHEP value

  // Sometimes, for interactions near threshold, Fermi momentum might bring
  // the neutrino energy in the nucleon rest frame below threshold (for the
  // selected interaction). In this case mark the event as unphysical and
  // abort the current thread.
  const KPhaseSpace & kps = interaction->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("FermiMover", pNOTICE)
                  << "Event below threshold after generating Fermi momentum";

     double Ethr = kps.Threshold();
     double Ev   = init_state->ProbeE(kRfHitNucRest);
     LOG("FermiMover", pNOTICE)
              << "Ev (@ nucleon rest frame) = " << Ev << ", Ethr = " << Ethr;

     evrec->EventFlags()->SetBitNumber(kBelowThrNRF, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E < Ethr after generating nucleon Fermi momentum");
     exception.SwitchOnFastForward();
     throw exception;
  }
  if(eject_nucleon_pdg != 0) {
    this->Emit2ndNucleonFromSRC(evrec, eject_nucleon_pdg);
  }
}
//___________________________________________________________________________
void FermiMover::Emit2ndNucleonFromSRC(GHepRecord * evrec,
                                       const int eject_pdg_code) const
{
  LOG("FermiMover", pINFO) << "Adding a recoil nucleon "
      "due to short range corelation";

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
  double M  = PDGLibrary::Instance()->Find(eject_pdg_code)->Mass();
  double E  = TMath::Sqrt(px*px+py*py+pz*pz+M*M);

  evrec->AddParticle(
      eject_pdg_code, status, imom, -1, -1, -1, px, py, pz, E, vx, vy, vz, 0);
}
//___________________________________________________________________________
void FermiMover::AddTargetNucleusRemnant(GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("FermiMover", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  for(int id = fd; id <= ld; id++) {

    // compute A,Z for final state nucleus & get its PDG code and its mass
    GHepParticle * particle = evrec->Particle(id);
    assert(particle);
    int  pdgc = particle->Pdg();
    bool is_p  = pdg::IsProton (pdgc);
    bool is_n  = pdg::IsNeutron(pdgc);

    if (is_p) Z--;
    if (is_p || is_n) A--;

    Px += particle->Px();
    Py += particle->Py();
    Pz += particle->Pz();
    E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn) {
    LOG("HadronicVtx", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("FermiMover", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);
}
//___________________________________________________________________________
void FermiMover::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FermiMover::LoadConfig(void)
{
  RgKey nuclkey = "NuclearModel";
  fNuclModel = 0;
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  this->GetParamDef("KeepHitNuclOnMassShell", fKeepNuclOnMassShell, false);
  this->GetParamDef("SimRecoilNucleon",       fSRCRecoilNucleon,    false);
}
//____________________________________________________________________________
