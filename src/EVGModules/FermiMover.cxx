//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 08, 2004

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
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "EVGModules/FermiMover.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/KineUtils.h"

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

  // check whether to emit a 2nd nucleon due to short range corelations
  if(fSRCRecoilNucleon) {
       this->Emit2ndNucleonFromSRC(evrec);
  }

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

  // generate a Fermi momentum & removal energy
  fNuclModel->GenerateNucleon(*tgt);
  TVector3 p3 = fNuclModel->Momentum3();
  double w    = fNuclModel->RemovalEnergy();

  LOG("FermiMover", pINFO) 
     << "Generated nucleon momentum: ("
     << p3.Px() << ", " << p3.Py() << ", " << p3.Pz() << "), "
     << "|p| = " << p3.Mag();
  LOG("FermiMover", pINFO) 
     << "Generated nucleon removal energy: w = " << w;
  
  double pF2 = p3.Mag2(); // (fermi momentum)^2

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  nucleon->SetRemovalEnergy(w);

  // struck nucleon energy:
  // two possible prescriptions depending on whether you want to force
  // the sruck nucleon to be on the mass-shell or not...

  double EN=0;

  if(!fKeepNuclOnMassShell) {
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
}
//___________________________________________________________________________
void FermiMover::Emit2ndNucleonFromSRC(GHepRecord * evrec) const
{
  LOG("FermiMover", pINFO) 
    << "Deciding whether to add a recoil nucleon "
    << "due to short range corelation";

  GHepParticle * nucleus = evrec->TargetNucleus();
  GHepParticle * nucleon = evrec->HitNucleon();

  // hit nuclear target & nucleon pdg codes
  int nucleus_pdgc = nucleus->Pdg();
  int nucleon_pdgc = nucleon->Pdg();

  // check the actual hit nucleon momentum momentum
  double pn = nucleon->P4()->Vect().Mag();

  // get kF
  string fKFTable = "Default";
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable(fKFTable);
  double kF = kft->FindClosestKF(nucleus_pdgc, nucleon_pdgc);

  // check whether to emit a nucleon due to short range corelation
  bool eject = (pn > kF);
  if(eject) {

     LOG("FermiMover", pDEBUG)
        << "pn = " << pn << " > kF = " << kF << " => " 
        << "Ejecting recoil nucleon";

    double Pp = (nucleon->Pdg() == kPdgProton) ? 0.05 : 0.95;

    RandomGen * rnd = RandomGen::Instance();
    double prob = rnd->RndGen().Rndm();
    int code = (prob > Pp) ? kPdgNeutron : kPdgProton;

    GHepStatus_t status = kIStHadronInTheNucleus;
    int imom = evrec->TargetNucleusPosition();

    //-- Has opposite momentum from the struck nucleon
    double vx = nucleon->Vx();
    double vy = nucleon->Vy();
    double vz = nucleon->Vz();
    double px = -1.* nucleon->Px();
    double py = -1.* nucleon->Py();
    double pz = -1.* nucleon->Pz();
    double M  = PDGLibrary::Instance()->Find(code)->Mass();
    double E  = TMath::Sqrt(px*px+py*py+pz*pz+M*M);

    evrec->AddParticle(
        code, status, imom, -1, -1, -1, px, py, pz, E, vx, vy, vz, 0);
  }
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
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fNuclModel = 0;

  RgKey nuclkey = "NuclearModel";
//  RgAlg nuclalg = fConfig->GetAlgDef(nuclkey, gc->GetAlg(nuclkey));
//  LOG("FermiMover", pINFO) << "Loading nuclear model: " << nuclalg;

  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);

  fKeepNuclOnMassShell = 
         fConfig->GetBoolDef("KeepHitNuclOnMassShell", false);

  fSRCRecoilNucleon = fConfig->GetBoolDef(
      "SimRecoilNucleon", gc->GetBool("SRC-SimRecoilNucleon"));
}
//____________________________________________________________________________

