//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGModules/QELHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
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
  //-- Determine the pdg & status code of the recoil nucleon
  Interaction * interaction = evrec->Summary();
  const Target & tgt = interaction->InitState().Tgt();
  int pdgc = interaction->RecoilNucleonPdg();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                           kIStHadronInTheNucleus : kIStStableFinalState;

  //-- Get nucleon 4-momentum (in the LAB frame) & position
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);
  TLorentzVector v4(0.,0.,0.,0.); 

  //-- Get mother position
  int mom = evrec->HitNucleonPosition();

  //-- Add the final state recoil nucleon at the EventRecord
  LOG("QELHadronicVtx", pINFO) << "Adding nucleon [pdgc = " << pdgc << "]";

  evrec->AddParticle(pdgc, ist, mom,-1,-1,-1, p4, v4);
}
//___________________________________________________________________________
