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
  //-- Determine the pdg code of the recoil nucleon
  Interaction * interaction = evrec->GetInteraction();
  int recoil_nuc_pdgc = utils::interaction::RecoilNucleonPdgCode(interaction);

  //-- Get all initial & final state particles 4-momenta (in the LAB frame)
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Add the final state recoil nucleon at the EventRecord
  LOG("QELHadronicVtx", pINFO)
                     << "Adding nucleon [pdgc = " << recoil_nuc_pdgc << "]";

  TLorentzVector v4(0.,0.,0.,0.);
  int mom = evrec->StruckNucleonPosition();

  evrec->AddParticle(recoil_nuc_pdgc,
                   kIStStableFinalState, mom,-1,-1,-1, p4, v4);
}
//___________________________________________________________________________
