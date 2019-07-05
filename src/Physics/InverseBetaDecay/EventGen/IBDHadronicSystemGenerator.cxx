//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl> - October 29, 2009
         using code from the QELKinematicGenerator written by
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 03, 2004

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/InverseBetaDecay/EventGen/IBDHadronicSystemGenerator.h"

using namespace genie;

//___________________________________________________________________________
IBDHadronicSystemGenerator::IBDHadronicSystemGenerator() :
HadronicSystemGenerator("genie::IBDHadronicSystemGenerator")
{

}
//___________________________________________________________________________
IBDHadronicSystemGenerator::IBDHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::IBDHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
IBDHadronicSystemGenerator::~IBDHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void IBDHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  // Add the recoil baryon 
  // (p or n)
  // Its 4-momentum is computed by requiring the energy + momentum to be
  // conserved.
  this->AddRecoilBaryon(evrec);
}
//___________________________________________________________________________
void IBDHadronicSystemGenerator::AddRecoilBaryon(GHepRecord * evrec) const
{
  //-- Determine the pdg & status code of the recoil baryon
  Interaction * interaction = evrec->Summary();
  const int pdgc = interaction->RecoilNucleonPdg();
  assert(pdgc!=0);

  //-- Determine the status code
  const Target & tgt = interaction->InitState().Tgt();
  const GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- Get the vtx position
  const GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  //-- Get nucleon 4-momentum (in the LAB frame) & position
  const TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Get mother position
  const int mom = evrec->HitNucleonPosition();

  //-- Add the final state recoil baryon at the EventRecord
  LOG("IBD", pINFO) 
      << "Adding recoil baryon [pdgc = " << pdgc << "]";

  GHepParticle p(pdgc, ist, mom,-1,-1,-1, p4, vtx);
  const double w = evrec->Particle(mom)->RemovalEnergy();
  p.SetRemovalEnergy(w);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
