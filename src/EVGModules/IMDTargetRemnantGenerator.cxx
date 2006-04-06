//____________________________________________________________________________
/*!

\class   genie::IMDTargetRemnantGenerator

\brief   Generates all the non - primarly lepton final state particles in v
         IMD interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 17, 2005

*/
//____________________________________________________________________________

#include "EVGModules/IMDTargetRemnantGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
IMDTargetRemnantGenerator::IMDTargetRemnantGenerator() :
EventRecordVisitorI("genie::IMDTargetRemnantGenerator")
{

}
//___________________________________________________________________________
IMDTargetRemnantGenerator::IMDTargetRemnantGenerator(string config) :
EventRecordVisitorI("genie::IMDTargetRemnantGenerator", config)
{

}
//___________________________________________________________________________
IMDTargetRemnantGenerator::~IMDTargetRemnantGenerator()
{

}
//___________________________________________________________________________
void IMDTargetRemnantGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this -> AddElectronNeutrino     (evrec);
  this -> AddTargetNucleusRemnant (evrec);
}
//___________________________________________________________________________
void IMDTargetRemnantGenerator::AddElectronNeutrino(GHepRecord * evrec) const
{
  //-- Get all initial & final state particles 4-momenta (in the LAB frame)

  //incoming v:
  GHepParticle * nu = evrec->Probe();

  //struck nucleon:
  GHepParticle * el = evrec->StruckElectron();

  //final state primary lepton:
  GHepParticle * l = evrec->FinalStatePrimaryLepton();

  assert(nu);
  assert(el);
  assert(l);

  //-- Force energy conservation

  // Pvmu(Ev,pxv,pyv,pzv) + Pe(En,pxn,pyn,pzn) = Pmu(El,pxl,pyl,pzl) + Pve
  double E  = nu->E()  + el->E()  - l->E();
  double px = nu->Px() + el->Px() - l->Px();
  double py = nu->Py() + el->Py() - l->Py();
  double pz = nu->Pz() + el->Pz() - l->Pz();

  //-- Add the final state recoil nucleon at the EventRecord

  LOG("IMDTargetRemnant", pINFO) << "Adding final state electron neutrino";

  int mom = evrec->StruckElectronPosition();
  evrec->AddParticle(
          kPdgNuE,kIStStableFinalState, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
void IMDTargetRemnantGenerator::AddTargetNucleusRemnant(
                                                    GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("IMDTargetRemnant", pDEBUG) << "Adding final state nucleus";

  //-- get A,Z for initial state nucleus
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  bool is_nucleus = init_state.GetTarget().IsNucleus();
  if (!is_nucleus) {
    LOG("IMDTargetRemnant", pDEBUG)
               << "Initial state not a nucleus - no remnant nucleus to add";
    return;
  }
  int A = init_state.GetTarget().A();
  int Z = init_state.GetTarget().Z();

  int    ipdgc = pdg::IonPdgCode(A, Z);
  double mass  = PDGLibrary::Instance()->Find(ipdgc)->Mass();

  //-- Add the nucleus to the event record
  LOG("IMDTargetRemnant", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                            << ", pdgc = " << ipdgc << "]";

  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
             ipdgc,kIStStableFinalState, mom,-1,-1,-1, 0,0,0,mass, 0,0,0,0);
}
//___________________________________________________________________________
