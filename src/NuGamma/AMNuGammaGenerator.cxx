//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 25, 2008 - CA
   This event generation modules was first added in version 2.3.1 as part of
   the new event generation thread handling amonaly-mediated single gamma
   interactions. 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "EVGModules/AMNuGammaGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;

//___________________________________________________________________________
AMNuGammaGenerator::AMNuGammaGenerator() :
EventRecordVisitorI("genie::AMNuGammaGenerator")
{

}
//___________________________________________________________________________
AMNuGammaGenerator::AMNuGammaGenerator(string config) :
EventRecordVisitorI("genie::AMNuGammaGenerator", config)
{

}
//___________________________________________________________________________
AMNuGammaGenerator::~AMNuGammaGenerator()
{

}
//___________________________________________________________________________
void AMNuGammaGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  this->AddPhoton(evrec);
  this->AddFinalStateNeutrino(evrec);
  this->AddTargetRemnant(evrec);
  this->AddRecoilNucleon(evrec);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddPhoton(GHepRecord * evrec) const
{
  LOG("AMNuGammaGenerator", pINFO) << "Adding final state photon";

  RandomGen * rnd = RandomGen::Instance();

  GHepParticle * nu = evrec->Probe();
  const TLorentzVector & vtx = *(nu->X4());
  const TLorentzVector & p4v = *(nu->P4());

  double Ev = p4v.Energy();
  double Eg = 0.5*Ev * (1 + rnd->RndKine().Rndm());
  TLorentzVector p4(p4v);
  p4 *= (Eg/Ev);

  GHepParticle p(kPdgGamma, kIStStableFinalState, 0,-1,-1,-1, p4, vtx);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const
{
  LOG("AMNuGammaGenerator", pINFO) << "Adding final state neutrino";

  GHepParticle * nu    = evrec->Probe(); // incoming v
  GHepParticle * gamma = evrec->Particle(nu->FirstDaughter()); // gamma
  assert(nu);
  assert(gamma);

  const TLorentzVector & p4nu    = *(nu->P4());
  const TLorentzVector & p4gamma = *(gamma ->P4());

  const TLorentzVector & vtx = *(nu->X4());
  TLorentzVector p4 = p4nu - p4gamma;
 
  GHepParticle p(nu->Pdg(), kIStStableFinalState, 0,-1,-1,-1, p4, vtx);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("AMNuGammaGenerator", pINFO) << "Adding final state nucleus";

  //-- skip for non nuclear targets
  GHepParticle * nucleus = evrec->TargetNucleus();
  if (!nucleus) {
    LOG("AMNuGammaGenerator", pDEBUG)
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
      LOG("AMNuGammaGenerator", pFATAL)
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
  LOG("AMNuGammaGenerator", pINFO)
       << "Adding nucleus [A = " << A << ", Z = " << Z
                                           << ", pdgc = " << ipdgc << "]";   
 
  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
           ipdgc,kIStStableFinalState, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddRecoilNucleon(GHepRecord * evrec) const
{
  LOG("AMNuGammaGenerator", pINFO) << "Adding recoil nucleon";

  //-- Get mother position
  int mom = evrec->HitNucleonPosition();

  //-- Determine the pdg & status code of the recoil baryon
  int pdgc = evrec->Particle(mom)->Pdg();
  assert(pdgc!=0);
  
  //-- Determine the status code
  Interaction * interaction = evrec->Summary();
  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                  kIStHadronInTheNucleus : kIStStableFinalState;
  
  //-- Get the vtx position
  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  //-- Get nucleon 4-momentum (in the LAB frame) & position
//  TLorentzVector p4 = this->Hadronic4pLAB(evrec);
  GHepParticle * hitnuc = evrec->HitNucleon();
  const TLorentzVector & p4n = *(hitnuc->P4());
  TLorentzVector p4(p4n);
  
 //-- Add the final state recoil baryon at the EventRecord
  LOG("AMNuGammaGenerator", pINFO)
      << "Adding recoil baryon [pdgc = " << pdgc << "]";

  GHepParticle p(pdgc, ist, mom,-1,-1,-1, p4, vtx);
  double w = evrec->Particle(mom)->RemovalEnergy();
  p.SetRemovalEnergy(w);
  evrec->AddParticle(p);
}
//___________________________________________________________________________


