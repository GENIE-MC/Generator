//____________________________________________________________________________
/*!

\class   genie::HadronicSystemGenerator

\brief   Abstract class. Is used to pass some commonly recurring methods to
         all concrete implementations of the EventRecordVisitorI interface
         generating the hadronic system for a specific processes (QEL,DIS,
         RES,...)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created Juy 16, 2005

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
// INtended for DIS hadronic system generators.

  TLorentzVector p4 = this->Hadronic4pLAB(evrec);
  TLorentzVector v4(0,0,0,0);

  int mom = evrec->StruckNucleonPosition();

  evrec->AddParticle(
        kPdgHadronicSyst, kIstDISPreFragmHadronicState, mom,-1,-1,-1, p4, v4);
}
//___________________________________________________________________________
void HadronicSystemGenerator::AddTargetNucleusRemnant(
                                                    GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("HadronicVtx", pDEBUG) << "Adding final state nucleus";

  //-- get A,Z for initial state nucleus

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  bool is_nucleus = init_state.GetTarget().IsNucleus();
  if (!is_nucleus) {
    LOG("HadronicVtx", pDEBUG)
               << "Initial state not a nucleus - no remnant nucleus to add";
    return;
  }
  int A = init_state.GetTarget().A();
  int Z = init_state.GetTarget().Z();

  //-- compute A,Z for final state nucleus & get its PDG code and its mass
  int  npdgc = init_state.GetTarget().StruckNucleonPDGCode();
  bool is_p  = pdg::IsProton(npdgc);
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
  double mass = particle->Mass();

  //-- Add the nucleus to the event record
  LOG("HadronicVtx", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                            << ", pdgc = " << ipdgc << "]";

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
             ipdgc,kIStStableFinalState, mom,-1,-1,-1, 0,0,0,mass, 0,0,0,0);
}
//___________________________________________________________________________
TLorentzVector HadronicSystemGenerator::Hadronic4pLAB(
                                                    GHepRecord * evrec) const
{
// Returns the final state hadronic system 4-p in LAB

  //incoming v:
  GHepParticle * nu = evrec->Probe();

  //struck nucleon:
  GHepParticle * N = evrec->StruckNucleon();

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

  double PX = nu->Px()     + N->Px()     - l->Px();
  double PY = nu->Py()     + N->Py()     - l->Py();
  double PZ = nu->Pz()     + N->Pz()     - l->Pz();
  double E  = nu->Energy() + N->Energy() - l->Energy();

  //TLorentzVector pX4 = (*nu->P4()) + (*N->P4()) - (*l->P4())

  TLorentzVector pX4(PX,PY,PZ,E);
  return pX4; 
}
//___________________________________________________________________________
TVector3 HadronicSystemGenerator::HCM2LAB(GHepRecord * evrec) const
{
// Velocity for the Hadronic CM -> LAB active Lorentz transform

  TLorentzVector pH = this->Hadronic4pLAB(evrec);

  //-- Compute the velocity of the LAB frame in the Final State Hadronic
  //   CM Frame (PxH/EH, PyH/EH, PzH/EH)

  TVector3 beta = pH.BoostVector();

  LOG("HadronicVtx", pINFO)
                   << "\n beta (HCM --> LAB): " << Vec3AsString(&beta);
  return beta;
}
//___________________________________________________________________________
int HadronicSystemGenerator::HadronShowerCharge(GHepRecord * evrec) const
{
// Returns the hadron shower charge in units of +e
// eg in v n -> l- X the hadron shower charge is +1

  int HadronShowerCharge = 0;

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  int hit_nucleon = init_state.GetTarget().StruckNucleonPDGCode();

  assert( pdg::IsProton(hit_nucleon) || pdg::IsNeutron(hit_nucleon) );

  double qfsl  = interaction->GetFSPrimaryLepton()->Charge() / 3.;
  double qinit = PDGLibrary::Instance()->Find(hit_nucleon)->Charge() / 3.;

  HadronShowerCharge = (int) (qinit - qfsl);

  return HadronShowerCharge;
}
//____________________________________________________________________________

