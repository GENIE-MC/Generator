//____________________________________________________________________________
/*!

\class   genie::QELHadronicSystemGenerator

\brief   Generates the final state hadronic system in v QEL interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include "EVGModules/QELHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepOrder.h"
#include "Interaction/IUtils.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
QELHadronicSystemGenerator::QELHadronicSystemGenerator() :
HadronicSystemGenerator("genie::QELHadronicSystemGenerator")
{

}
//___________________________________________________________________________
QELHadronicSystemGenerator::QELHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::QELHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
QELHadronicSystemGenerator::~QELHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void QELHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- If the struck nucleon was within a nucleus, then add the final state
  //   nucleus at the EventRecord
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the recoil nucleon
  //   Its 4-momentum is computed by requiring the energy + momentum to be
  //   conserved.
  this->AddRecoilNucleon(evrec);
}
//___________________________________________________________________________
void QELHadronicSystemGenerator::AddRecoilNucleon(GHepRecord * evrec) const
{
  //-- Get the interaction & initial state objects

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  //-- Determine the pdg code of the recoil nucleon

  int recoil_nuc_pdgc = utils::interaction::RecoilNucleonPdgCode(interaction);

  //-- Get all initial & final state particles 4-momenta (in the LAB frame)

  //incoming v:
  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfLab);
  assert(nu_p4);

  //struck nucleon:
  TLorentzVector * nucl_p4 = init_state.GetTarget().StruckNucleonP4();
  assert(nucl_p4);

  //final state primary lepton:
  int fsl_pdgc = interaction->GetFSPrimaryLepton()->PdgCode();
  GHepParticle * fsl = evrec->FindParticle(fsl_pdgc, kIStStableFinalState, 0);
  assert(fsl);

  //-- Force energy conservation

  // Pv(Ev,pxv,pyv,pzv) + Pnucl(En,pxn,pyn,pzn) = Pl(El,pxl,pyl,pzl) + Precoil

  double E  = nu_p4->Energy() + nucl_p4->Energy() - fsl->E();
  double px = nu_p4->Px()     + nucl_p4->Px()     - fsl->Px();
  double py = nu_p4->Py()     + nucl_p4->Py()     - fsl->Py();
  double pz = nu_p4->Pz()     + nucl_p4->Pz()     - fsl->Pz();

  delete nu_p4;

  //-- Add the final state recoil nucleon at the EventRecord

  LOG("QELHadronicVtx", pINFO)
                     << "Adding nucleon [pdgc = " << recoil_nuc_pdgc << "]";

  int mom = GHepOrder::StruckNucleonPosition(interaction);

  evrec->AddParticle(recoil_nuc_pdgc,
                   kIStStableFinalState, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
