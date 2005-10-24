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
#include "GHEP/GHepOrder.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

using namespace genie;

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
  double mass  = particle->Mass();

  //-- Add the nucleus to the event record
  LOG("HadronicVtx", pINFO)
          << "Adding nucleus [A = " << A << ", Z = " << Z
                                            << ", pdgc = " << ipdgc << "]";

  int mom = GHepOrder::TargetNucleusPosition(interaction);

  evrec->AddParticle(
             ipdgc,kIStStableFinalState, mom,-1,-1,-1, 0,0,0,mass, 0,0,0,0);
}
//___________________________________________________________________________
