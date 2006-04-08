//____________________________________________________________________________
/*!

\class   genie::COHHadronicSystemGenerator

\brief   Generates the final state hadronic system in v COH interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <cstdlib>

#include "EVGModules/COHHadronicSystemGenerator.h"
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
COHHadronicSystemGenerator::COHHadronicSystemGenerator() :
HadronicSystemGenerator("genie::COHHadronicSystemGenerator")
{

}
//___________________________________________________________________________
COHHadronicSystemGenerator::COHHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::COHHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
COHHadronicSystemGenerator::~COHHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void COHHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  //-- Add the final state nucleus 
  this->AddTargetNucleusRemnant(evrec);

  //-- Add the final state pion
  this->AddFinalStatePion(evrec);
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::AddTargetNucleusRemnant(
                                                   GHepRecord * evrec) const
{
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleus);

  int mom   = evrec->TargetNucleusPosition();
  int ipdgc = nucleus->PdgCode();

  evrec->AddParticle(
     ipdgc,kIStStableFinalState, mom,-1,-1,-1, 
         nucleus->Px(), nucleus->Py(), nucleus->Pz(), nucleus->E(), 0,0,0,0);
}
//___________________________________________________________________________
void COHHadronicSystemGenerator::AddFinalStatePion(GHepRecord * evrec) const
{
  //-- Determine the pdg code of the final state pion
  Interaction * interaction = evrec->GetInteraction();
  const XclsTag & xcls_tag  = interaction->GetExclusiveTag();

  int pion_pdgc = 0;
  if      (xcls_tag.NPi0()     ==1) pion_pdgc = kPdgPi0;
  else if (xcls_tag.NPiPlus()  ==1) pion_pdgc = kPdgPiPlus;
  else if (xcls_tag.NPiMinus() ==1) pion_pdgc = kPdgPiMinus;
  else {
     LOG("COHHadronicVtx", pFATAL)
                          << "No final state pion information in XclsTag!";
     exit(1);
  }

  //-- Figure out the final state pion 4-p
  TLorentzVector p4(0,0,0,0);

  //-- Figure out the mother particle position
  int mom = evrec->TargetNucleusPosition();
  assert(mom>0);

  //-- Add the final state pion at the EventRecord
  LOG("COHHadronicVtx", pINFO)
                             << "Adding pion [pdgc = " << pion_pdgc << "]";
  TLorentzVector v4(0.,0.,0.,0.);

  evrec->AddParticle(pion_pdgc, kIStStableFinalState, mom,-1,-1,-1, p4, v4);
}
//___________________________________________________________________________

