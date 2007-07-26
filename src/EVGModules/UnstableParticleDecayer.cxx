//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - November 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <algorithm>
#include <sstream>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TParticlePDG.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Decay/DecayModelI.h"
#include "EVGModules/UnstableParticleDecayer.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using std::count;
using std::ostringstream;

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer() :
EventRecordVisitorI("genie::UnstableParticleDecayer")
{
  fDecayers = 0;
}
//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer(string config) :
EventRecordVisitorI("genie::UnstableParticleDecayer", config)
{
  fDecayers = 0;
}
//___________________________________________________________________________
UnstableParticleDecayer::~UnstableParticleDecayer()
{
  if(fDecayers) {
    fDecayers->clear();
    delete fDecayers;
  }
}
//___________________________________________________________________________
void UnstableParticleDecayer::ProcessEventRecord(GHepRecord * evrec) const
{
// Loop over particles, find unstable ones and decay them

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  unsigned int ipos = 0;

  vector<const DecayModelI *>::const_iterator dec_iter;

  //-- Check whether the interaction is off a nuclear target or free nucleon
  //   Depending on whether this module is run before or after the hadron
  //   transport module it would affect the daughters status code
  GHepParticle * nucltgt = evrec->TargetNucleus();
  bool in_nucleus = (nucltgt!=0);

  while( (p = (GHepParticle *) piter.Next()) ) {

     if( this->ToBeDecayed(p) ) {
        LOG("ParticleDecayer", pNOTICE)
                   << "Decaying unstable particle: " << p->Name();

        //-- find the first decayer to handle the current particle
        dec_iter = fDecayers->begin();
        fCurrDecayer=0;
        for( ; dec_iter != fDecayers->end(); ++dec_iter) {
          const DecayModelI * decayer = *dec_iter;
          LOG("ParticleDecayer", pNOTICE)
                 << "Requesting decay from " << decayer->Id().Key();
          if(decayer->IsHandled(p->Pdg())) {
            fCurrDecayer = decayer;
            break;
          }
        }
        //-- handle the case where no decayer is found
        if(fCurrDecayer==0) {
           LOG("ParticleDecayer", pWARN) 
            << "Couldn't find a decayer for: " << p->Name() << ". Skipping!";
           continue;        
        }

        //-- Decay it & retrieve the decay products
        //   The decayer might not be able to handle it - in which case it
        //   should return a NULL decay products container

        TLorentzVector p4 = *(p->P4());
        DecayerInputs_t dinp;
        dinp.PdgCode = p->Pdg();
        dinp.P4      = &p4;

        TClonesArray * decay_products = fCurrDecayer->Decay(dinp);

        //-- Check whether the particle was decayed
        if(decay_products) {
           LOG("ParticleDecayer", pINFO) << "The particle was decayed";

           //-- Mark it as a 'decayed state' & add its daughter links
           p->SetStatus(kIStDecayedState);

           //-- Add the mom & daughters to the event record
           this->CopyToEventRecord(decay_products, evrec, p, ipos, in_nucleus);

           //-- Update the event weight for each weighted particle decay
           double decay_weight = fCurrDecayer->Weight();
           evrec->SetWeight(evrec->Weight() * decay_weight);

           //-- Clean-up decay products
           decay_products->Delete();
           delete decay_products;
        }// !=0
     }// to be decayed?
     ipos++;

  } // loop over particles

  LOG("ParticleDecayer", pNOTICE) 
          << "Done finding unstable particles & decaying them!";
}
//___________________________________________________________________________
bool UnstableParticleDecayer::ToBeDecayed(GHepParticle * particle) const
{
   GHepStatus_t ist = (fRunBefHadroTransp) ?
                       kIStHadronInTheNucleus : kIStStableFinalState;

   if( particle->Pdg() != 0 && particle->Status() == ist)
                                   return this->IsUnstable(particle);
   return false;
}
//___________________________________________________________________________
bool UnstableParticleDecayer::IsUnstable(GHepParticle * particle) const
{
  int pdg_code = particle->Pdg();

  TParticlePDG * ppdg = PDGLibrary::Instance()->Find(pdg_code);
  if( ppdg->Lifetime() < fMaxLifetime ) { /* ... */ };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - temporary/
  // ROOT's TParticlepdg::Lifetime() does not work properly
  // do something else instead (temporarily)

  if( fRunBefHadroTransp ) {
    /* run before the hadron transport MC */
    if(utils::res::IsBaryonResonance(pdg_code)) return true;

  } else {
    /* run after the hadron transport MC - decay particles with ct<~0.1mm*/   
    int particles_to_decay[] = {
          kPdgPi0, 
          kPdgEta, kPdgEtaPrm, 
          kPdgRho0, kPdgRhoP, kPdgRhoM,
          kPdgomega,kPdgPhi };

    const int N = sizeof(particles_to_decay) / sizeof(int);
    int matches = count(particles_to_decay, particles_to_decay+N, pdg_code);

    if(matches > 0 || utils::res::IsBaryonResonance(pdg_code)) return true;
  } 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - /temporary

  return false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::CopyToEventRecord(
  TClonesArray * decay_products, GHepRecord * evrec, GHepParticle * p,
                                       int mother_pos, bool in_nucleus) const
{
// Adds the decayer output at the GHEP event record. 
// Only adds the decay products - the mother particle already exists and is
// currently not duplicated. The mother's daughter list will be automatically 
// updated with every AddParticle() call

  LOG("ParticleDecayer", pINFO) << "Copying decay to event record...";

  int new_mother_pos = mother_pos;
  TMCParticle * dpmc =  0;

  TLorentzVector parent_x4 = *(p->X4());

  TObjArrayIter decay_iter(decay_products);

  while( (dpmc = (TMCParticle *) decay_iter.Next()) ) {

     int pdg = dpmc->GetKF();
     GHepStatus_t ist = GHepStatus_t (dpmc->GetKS()); 

/*
     double px = dpmc->GetPx();
     double py = dpmc->GetPy();
     double pz = dpmc->GetPz();
     double E  = dpmc->GetEnergy();
*/
     TLorentzVector p4(dpmc->GetPx(), dpmc->GetPy(), dpmc->GetPz(), dpmc->GetEnergy()); 

     TLorentzVector x4(dpmc->GetVx(), dpmc->GetVy(), dpmc->GetVz(), dpmc->GetTime()); 
     x4 += parent_x4;

     //-- and now add the decay products
     if(ist == kIStStableFinalState) { 

        LOG("ParticleDecayer", pDEBUG) << "Adding daughter... PDG=" << pdg;

        // figure out the right status code for the current daughter
        bool is_hadron = pdg::IsHadron(pdg);
        bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

        GHepStatus_t istfin = (hadron_in_nuc) ?
                               kIStHadronInTheNucleus : kIStStableFinalState;

        evrec->AddParticle(pdg, istfin, new_mother_pos,-1,-1,-1, p4, x4);
     }
  }
}
//___________________________________________________________________________
void UnstableParticleDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void UnstableParticleDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void UnstableParticleDecayer::LoadConfig(void)
{
  //-- Get the specified maximum lifetime tmax (decay with lifetime < tmax)
  //
  fMaxLifetime = fConfig->GetDoubleDef("MaxLifetime", 1e-9);

  //-- Check whether the module is being run before or after the hadron
  //   transport (intranuclear rescattering) module.
  //
  //   If it is run before the hadron transport (and after the hadronization)
  //   step it should decay only "unstable" particles (marked as hadrons in
  //   the nucleus) which would typically decay within the time required to 
  //   exit the nucleus - so, the algorithm wouldn't decay particles that 
  //   have to be rescattered first. In case that the generated event is off
  //   a free nucleon target, thi instance of the algorithm should do nothing.
  //
  //   If it is run after the hadon transport, then it should decay all the 
  //   'unstable' particles marked as 'present in the final state' and which 
  //   should be decay before the event is passed to the detector particle
  //   transport MC.
  //
  fRunBefHadroTransp = fConfig->GetBool("RunBeforeHadronTransport");

  //-- Load particle decayers
  //   Order is important if both decayers can handle a specific particle
  //   as only the first would get the chance to decay it

  int ndec = fConfig->GetIntDef("NDecayers",0);
  assert(ndec>0);

  if(fDecayers) {
    fDecayers->clear();
    delete fDecayers;
  }
  fDecayers = new vector<const DecayModelI *>(ndec);

  for(int idec = 0; idec < ndec; idec++) {
     ostringstream alg_key;
     alg_key     << "Decayer-" << idec;
     const DecayModelI * decayer = 
              dynamic_cast<const DecayModelI *> (this->SubAlg(alg_key.str()));
     (*fDecayers)[idec] = decayer;
  }
}
//___________________________________________________________________________
