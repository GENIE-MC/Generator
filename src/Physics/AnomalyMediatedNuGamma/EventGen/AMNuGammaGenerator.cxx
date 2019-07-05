//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 25, 2008 - CA
   This event generation modules was first added in version 2.3.1 as part of
   the new event generation thread handling amonaly-mediated single gamma
   interactions. 
 @ Sep 25, 2008 - CA
   Improved the calculation of the photon 4-momentum. Generating the 4-momentum
   at NRF' and then rotating it to NRF and boosting it back at the LAB.
   The photon cos(theta) follows a uniform distribution (with respect to the
   incoming neutrnino in NRF). Need further inputs for modeling the energy 
   transfer to the photon (uniform distribution up to the available energy).
   The final state neutrino 4-momentum determined from energy conservation
   (incoming nu = outgoing gamma + outgoing nu) but difficult to keep on the
   mass shell. The nucleon recoil is negligible. For nuclear targets the hit
   nucleon is forced back on the mass shell before intranuclear rescattering.
   The necessary energy is taken from the remnant nucleus.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/AnomalyMediatedNuGamma/EventGen/AMNuGammaGenerator.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;
using namespace genie::constants;

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
  this->AddRecoilNucleon(evrec);
//this->AddTargetRemnant(evrec);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddPhoton(GHepRecord * evrec) const
{
// Adding the final state photon
//
  LOG("AMNuGammaGenerator", pINFO) << "Adding final state photon";

  RandomGen * rnd = RandomGen::Instance();

  // Get boost vector for transforms between LAB <-> NRF (nucleon rest frame)
  GHepParticle * nuc = evrec->HitNucleon();
  const TLorentzVector & p4nuc_lab = *(nuc->P4());
  TVector3 beta = p4nuc_lab.BoostVector();

  // Get the neutrino 4-momentum at the LAB
  GHepParticle * nu = evrec->Probe();
  const TLorentzVector & p4v_lab = *(nu->P4()); 

  // Get the neutrino 4-momentum at the NRF
  TLorentzVector p4v_nrf = p4v_lab; 
  p4v_nrf.Boost(-1.*beta);                  

  // Generate the photon cos(theta) with respect to the neutrino direction
  // (at NRF) and a uniform azimuthal angle phi
  double costheta_gamma = -1.0 + 2.0 * rnd->RndKine().Rndm();
  double phi_gamma      =  2.0 * kPi * rnd->RndKine().Rndm();

  // Generate the fraction of the neutrino energy taken by the photon
  double efrac_gamma =  rnd->RndKine().Rndm();

  // Calculate the photon energy at the NRF
  double Ev_nrf = p4v_nrf.Energy();
  double Egamma_nrf = Ev_nrf * efrac_gamma;

  // Calculate the photon momentum components at a rotated NRF (NRF') 
  // where z is along the neutrino direction 
  double sintheta_gamma = TMath::Sqrt(1-TMath::Power(costheta_gamma,2));
  double pgamma_nrf_p = Egamma_nrf * costheta_gamma; // p(//)
  double pgamma_nrf_t = Egamma_nrf * sintheta_gamma; // p(-|)
  double px = pgamma_nrf_t * TMath::Sin(phi_gamma);
  double py = pgamma_nrf_t * TMath::Cos(phi_gamma);
  double pz = pgamma_nrf_p;

  // Take a unit vector along the neutrino direction @ the NRF
  TVector3 unit_nudir = p4v_nrf.Vect().Unit(); 

  // Rotate the photon momentum vector from NRF' -> NRF
  TVector3 p3gamma_nrf(px,py,pz);
  p3gamma_nrf.RotateUz(unit_nudir);

  // Get the photon 4-momentum back at the LAB
  TLorentzVector p4gamma_lab(p3gamma_nrf, Egamma_nrf);
  p4gamma_lab.Boost(beta); 

  // Add the photon at the event record
  const TLorentzVector & vtx = *(nu->X4());
  GHepParticle p(kPdgGamma,kIStStableFinalState,0,-1,-1,-1,p4gamma_lab,vtx);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const
{
// Adding the final state neutrino
// Just use 4-momentum conservation (init_neutrino = photon + final_neutrino)

  LOG("AMNuGammaGenerator", pINFO) << "Adding final state neutrino";

  GHepParticle * nu    = evrec->Probe();                       // incoming v
  GHepParticle * gamma = evrec->Particle(nu->FirstDaughter()); // gamma
  assert(nu);
  assert(gamma);

  const TLorentzVector & vtx = *(nu->X4()); // vtx

  const TLorentzVector & p4nu_lab    = *(nu->P4());
  const TLorentzVector & p4gamma_lab = *(gamma->P4());
  TLorentzVector p4_lab = p4nu_lab - p4gamma_lab;
 
  GHepParticle p(nu->Pdg(), kIStStableFinalState, 0,-1,-1,-1, p4_lab, vtx);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddRecoilNucleon(GHepRecord * evrec) const
{
// Adding the recoil nucleon.
// Recoil is negligible at the H3 model. However, for nuclear targets, the 
// hit nucleon was off-the-mass-shell (bound). Doing some magic here to bring 
// it back on-the-mass-shell so that it can be INTRANUKE'ed and appear in 
// the final state. The nucleon will keep its original (Fermi) 3-momentum but 
// take some energy back from the remnant nucleus. No such tweaking takes
// place for free nucleon targets.

  LOG("AMNuGammaGenerator", pINFO) << "Adding recoil nucleon";

  // Get the interaction target
  GHepParticle * tgt_nucleus = evrec->TargetNucleus();
  bool is_nuclear_target = (tgt_nucleus!=0);

  // Get the hit nucleon
  GHepParticle * hitnuc = evrec->HitNucleon();
  assert(hitnuc);
  
  // Get the hit nucleon pdg code (= recoil nucleon pdg code)
  int pdgc = hitnuc->Pdg();

  // Get the hit nucleon 4-momentum (LAB) 
  const TLorentzVector & p4n = *(hitnuc->P4());
  TLorentzVector p4(p4n);

  // Tweak the 4-momentum to bring the recoil nucleon on-the-mass-shell
  // (for nuclear targets only)
  if (is_nuclear_target) {
    double p = p4.Vect().Mag();
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    double E = TMath::Sqrt(m*m+p*p);
    p4.SetE(E);
  }

  // Get the vtx position
  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  // Add the recoil nucleon at the event record
  GHepStatus_t ist = (is_nuclear_target) ?
                      kIStHadronInTheNucleus : kIStStableFinalState;
  int mom = evrec->HitNucleonPosition();

  LOG("AMNuGammaGenerator", pINFO)
                  << "Adding recoil baryon [pdgc = " << pdgc << "]";

  GHepParticle p(pdgc, ist, mom,-1,-1,-1, p4, vtx);
//  double w = hitnuc->RemovalEnergy();
///  p.SetRemovalEnergy(w);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void AMNuGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
// Add the remnant nuclear target at the GHEP record

  LOG("AMNuGammaGenerator", pINFO) << "Adding final state nucleus";

  // Skip for non nuclear targets
  GHepParticle * nucleus = evrec->TargetNucleus();
  if (!nucleus) {
    LOG("AMNuGammaGenerator", pDEBUG)
     << "No nucleus in the initial state - no remnant nucleus to add in the f/s";
    return;
  }

  // Compute A,Z for final state nucleus & get its PDG code and its mass
  //GHepParticle * nucleon = evrec->HitNucleon();

  GHepParticle * hit_nucleon = evrec->HitNucleon(); // hit
  GHepParticle * rec_nucleon = evrec->Particle(hit_nucleon->FirstDaughter()); // recoil

  assert(rec_nucleon);
  int  npdgc = rec_nucleon->Pdg();
  bool is_p  = pdg::IsProton(npdgc);
  int A = nucleus->A();
  int Z = nucleus->Z();
  if (is_p) Z--;
  A--;
  int ipdgc = pdg::IonPdgCode(A, Z);
  TParticlePDG * remnant = PDGLibrary::Instance()->Find(ipdgc);
  if(!remnant) {
      LOG("AMNuGammaGenerator", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
      assert(remnant);
  }

  // Figure out the remnant 4-momentum
  double Mf  = remnant->Mass();   
  double Mf2 = TMath::Power(Mf,2);
  double px  = -1.* rec_nucleon->Px();  
  double py  = -1.* rec_nucleon->Py();
  double pz  = -1.* rec_nucleon->Pz();
  double E   = TMath::Sqrt(Mf2 + rec_nucleon->P4()->Vect().Mag2());
  E += (hit_nucleon->P4()->E() - rec_nucleon->P4()->E());
   
  //-- Add the nucleus to the event record
  LOG("AMNuGammaGenerator", pINFO)
       << "Adding nucleus [A = " << A << ", Z = " << Z
                                           << ", pdgc = " << ipdgc << "]";   
  int mom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
           ipdgc,kIStStableFinalState, mom,-1,-1,-1, px,py,pz,E, 0,0,0,0);
}
//___________________________________________________________________________
