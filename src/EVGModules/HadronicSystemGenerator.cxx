//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 16, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/HadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::utils::print;

//___________________________________________________________________________
HadronicSystemGenerator::HadronicSystemGenerator() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
HadronicSystemGenerator::HadronicSystemGenerator(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
HadronicSystemGenerator::HadronicSystemGenerator(string name, string config):
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
HadronicSystemGenerator::~HadronicSystemGenerator()
{

}
//___________________________________________________________________________
void HadronicSystemGenerator::AddFinalHadronicSyst(GHepRecord * evrec) const
{
// Adds a GHEP entry for the sum of the f/s hadronic system.
// Intended for DIS hadronic system generators.

  TLorentzVector p4 = this->Hadronic4pLAB(evrec);
  TLorentzVector v4(0,0,0,0);

  int mom = evrec->HitNucleonPosition();

  evrec->AddParticle(
       kPdgHadronicSyst, kIStDISPreFragmHadronicState, mom,-1,-1,-1, p4, v4);

  // update the interaction summary
  evrec->Summary()->KinePtr()->SetHadSystP4(p4);
}
//___________________________________________________________________________
void HadronicSystemGenerator::AddTargetNucleusRemnant(
                                                    GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("HadronicVtx", pDEBUG) << "Adding final state nucleus";

  //-- skip for non nuclear targets

  GHepParticle * nucleus = evrec->TargetNucleus();
  if (!nucleus) {
    LOG("HadronicVtx", pDEBUG)
               << "Initial state not a nucleus - no remnant nucleus to add";
    return;
  }

  //-- compute A,Z for final state nucleus & get its PDG code and its mass

  GHepParticle * nucleon = evrec->HitNucleon();
  assert(nucleon);
  int  npdgc = nucleon->Pdg();
  bool is_p  = pdg::IsProton(npdgc);
  int A = nucleus->A();
  int Z = nucleus->Z();
  if (is_p) Z--;
  A--;

  TParticlePDG * particle = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  particle = PDGLibrary::Instance()->Find(ipdgc);
  if(!particle) {
      LOG("HadronicVtx", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
      assert(particle);
  }
  double Mf  = particle->Mass();   // remnant nucleus rest mass
  double Mf2 = TMath::Power(Mf,2);

  //-- Has opposite momentum from the struck nucleon

  double px = -1.* nucleon->Px();
  double py = -1.* nucleon->Py();
  double pz = -1.* nucleon->Pz();
  double E  = TMath::Sqrt(Mf2 + nucleon->P4()->Vect().Mag2());

  //-- Add the nucleus to the event record
  LOG("HadronicVtx", pINFO)
       << "Adding nucleus [A = " << A << ", Z = " << Z
                                           << ", pdgc = " << ipdgc << "]";

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
           ipdgc,kIStStableFinalState, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
TLorentzVector HadronicSystemGenerator::Hadronic4pLAB(
                                                    GHepRecord * evrec) const
{
// Returns the final state hadronic system 4-p in LAB

  //incoming v:
  GHepParticle * nu = evrec->Probe();

  //struck nucleon:
  GHepParticle * N = evrec->HitNucleon();

  //final state primary lepton:
  GHepParticle * l = evrec->FinalStatePrimaryLepton();

  assert(nu);
  assert(N);
  assert(l);

  LOG("HadronicVtx", pINFO)
      << "\n v [LAB]: " << P4AsString( nu->P4() )
      << "\n N [LAB]: " << P4AsString( N->P4()  )
      << "\n l [LAB]: " << P4AsString( l->P4()  );

  //-- Compute the Final State Hadronic System 4p (PX = Pv + PN - Pl)

  //double PX = nu->Px()     + N->Px()     - l->Px();
  //double PY = nu->Py()     + N->Py()     - l->Py();
  //double PZ = nu->Pz()     + N->Pz()     - l->Pz();
  //double E  = nu->Energy() + N->Energy() - l->Energy();
  //TLorentzVector pX4(PX,PY,PZ,E);

  const TLorentzVector & p4nu = *(nu->P4());
  const TLorentzVector & p4N  = *(N ->P4());
  const TLorentzVector & p4l  = *(l ->P4());

  TLorentzVector pX4 = p4nu + p4N - p4l;

  LOG("HadronicVtx", pNOTICE) << "\n HadrSyst [LAB]: " << P4AsString(&pX4);

  return pX4; 
}
//___________________________________________________________________________
TLorentzVector HadronicSystemGenerator::MomentumTransferLAB(
                                                    GHepRecord * evrec) const
{
  //incoming v:
  GHepParticle * nu = evrec->Probe();

  //final state primary lepton:
  GHepParticle * l = evrec->FinalStatePrimaryLepton();

  assert(nu);
  assert(l);

  const TLorentzVector & p4nu = *(nu->P4());
  const TLorentzVector & p4l  = *(l ->P4());

  TLorentzVector pq4 = p4nu - p4l; // q

  LOG("HadronicVtx", pNOTICE) 
                      << "\n Momentum Transfer [LAB]: " << P4AsString(&pq4);
  return pq4; 
}
//___________________________________________________________________________
TVector3 HadronicSystemGenerator::HCM2LAB(GHepRecord * evrec) const
{
// Velocity for the Hadronic CM -> LAB active Lorentz transform

  TLorentzVector pH = this->Hadronic4pLAB(evrec);

  //-- Compute the velocity of the LAB frame in the Final State Hadronic
  //   CM Frame (PxH/EH, PyH/EH, PzH/EH)

  TVector3 beta = pH.BoostVector();

  LOG("HadronicVtx", pINFO) << "beta (HCM->LAB): " << Vec3AsString(&beta);

  return beta;
}
//___________________________________________________________________________
int HadronicSystemGenerator::HadronShowerCharge(GHepRecord * evrec) const
{
// Returns the hadron shower charge in units of +e
// eg in v n -> l- X the hadron shower charge is +1

  int HadronShowerCharge = 0;

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  int hit_nucleon = init_state.Tgt().HitNucPdg();

  assert( pdg::IsProton(hit_nucleon) || pdg::IsNeutron(hit_nucleon) );

  double qfsl  = interaction->FSPrimLepton()->Charge() / 3.;
  double qinit = PDGLibrary::Instance()->Find(hit_nucleon)->Charge() / 3.;

  HadronShowerCharge = (int) (qinit - qfsl);

  return HadronShowerCharge;
}
//____________________________________________________________________________

