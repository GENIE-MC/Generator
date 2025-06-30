//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/EventGen/QELHadronicSystemGenerator.h"

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

  // Add the recoil baryon
  // (p or n - Lambda_c+,Sigma_c+,Sigma_c++ in charm/QEL)
  // Its 4-momentum is computed by requiring the energy + momentum to be
  // conserved.
  this->AddRecoilBaryon(evrec);
}
//___________________________________________________________________________
void QELHadronicSystemGenerator::AddRecoilBaryon(GHepRecord * evrec) const
{
  //-- Determine the pdg & status code of the recoil baryon
  Interaction * interaction = evrec->Summary();
  const XclsTag & xcls = interaction->ExclTag();
  int pdgc = 0;
  if     (xcls.IsCharmEvent())   { pdgc = xcls.CharmHadronPdg();           }
  else if(xcls.IsStrangeEvent()) { pdgc = xcls.StrangeHadronPdg();         }
  else                           { pdgc = interaction->RecoilNucleonPdg(); }
  assert(pdgc!=0);

  //-- Determine the status code
  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- Get the vtx position
  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  //-- Get nucleon 4-momentum (in the LAB frame) & position
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Get mother position
  int mom = evrec->HitNucleonPosition();

  //-- Add the final state recoil baryon at the EventRecord
  LOG("QELHadronicVtx", pINFO)
      << "Adding recoil baryon [pdgc = " << pdgc << "]";

  GHepParticle p(pdgc, ist, mom,-1,-1,-1, p4, vtx);
  double w = ( xcls.IsCharmEvent() || xcls.IsStrangeEvent()) ?
                  0. : evrec->Particle(mom)->RemovalEnergy();

  p.SetRemovalEnergy(w);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
