//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 05, 2008 - CA
   Added option to force charmed hadron decays
 @ Oct 10, 2008 - CA
   Added option to inhibit pi0 decays. Convert all particle positions to fm
   before adding them at the event record.
 @ Sep 30, 2009 - CA
   Requested to allow tau decays. Instead of providing yet another flag, 
   I am now allowing the user to specify a complete list of particles to
   be decayed by specifying a series of 
   <param type="bool" name="DecayParticleWithCode=[pdgcode]"> true </param>  
   parameters at UserPhysicsOptions.xml 
 @ Oct 02, 2009 - CA
   Requested to allow inhibiting decay channels. This is now possible by
   specifying a series of 
   <param type="bool" name="InhibitDecay/Particle=[code],Channel=[code]"> 
   true </param> options at UserPhysicsOptions.xml.
   Note that enabling this option results in **weighted** events.
   Solved problem with, say, inhibiting pi0 decay at this module, but having 
   other pi0 decayed deep in pythia when it decays hadrons (eg rho0) having
   pi0 in their decay products.

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

#include "Algorithm/AlgConfigPool.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Decay/DecayModelI.h"
#include "EVGModules/UnstableParticleDecayer.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Registry/Registry.h"
#include "Utils/StringUtils.h"

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

     LOG("ParticleDecayer", pINFO) << "Checking: " << p->Name();

     if( this->ToBeDecayed(p) ) {
        LOG("ParticleDecayer", pINFO)
           << "Decaying unstable particle: " << p->Name();

        //-- find the first decayer to handle the current particle
        dec_iter = fDecayers->begin();
        fCurrDecayer=0;
        for( ; dec_iter != fDecayers->end(); ++dec_iter) {
          const DecayModelI * decayer = *dec_iter;
          LOG("ParticleDecayer", pINFO)
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

  LOG("ParticleDecayer", pINFO) 
          << "Done finding unstable particles & decaying them!";
}
//___________________________________________________________________________
bool UnstableParticleDecayer::ToBeDecayed(GHepParticle * particle) const
{
   if(particle->Pdg() != 0) {
     bool check = false;
     GHepStatus_t ist = particle->Status();

     if(fRunBefHadroTransp) {
       check = (ist == kIStHadronInTheNucleus || ist == kIStStableFinalState);
     } else {
       check = (ist == kIStStableFinalState);
     }
     if(check) { return this->IsUnstable(particle); }
   }
   return false;
}
//___________________________________________________________________________
bool UnstableParticleDecayer::IsUnstable(GHepParticle * particle) const
{
  int pdg_code = particle->Pdg();

  // ROOT's TParticlepdg::Lifetime() does not work properly
  // do something else instead (temporarily)
  //
  // TParticlePDG * ppdg = PDGLibrary::Instance()->Find(pdg_code);
  //if( ppdg->Lifetime() < fMaxLifetime ) { /* ... */ };
  //

  // <temp/>
  if( fRunBefHadroTransp ) {
    //
    // Run *before* the hadron transport MC 
    // At this point we decay only baryon resonances
    //
    bool decay = utils::res::IsBaryonResonance(pdg_code);
    return decay;
  } 
  else {
    //
    // Run *after* the hadron transport MC 
    // At this point we decay only particles in the fParticlesToDecay 
    // PDGCodeList (filled in from config inputs)
    //
    bool decay = fParticlesToDecay.ExistsInPDGCodeList(pdg_code);
    return decay;
  } 
  // </temp>

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

     TLorentzVector p4(dpmc->GetPx(), 
                       dpmc->GetPy(), 
                       dpmc->GetPz(), 
                       dpmc->GetEnergy()); 
     TLorentzVector x4(dpmc->GetVx() / units::fm, 
                       dpmc->GetVy() / units::fm, 
                       dpmc->GetVz() / units::fm, 
                       0); 
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

  fAllowReconfig = false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Load particle decayers
  // Order is important if both decayers can handle a specific particle
  // as only the first would get the chance to decay it

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

  // Get the specified maximum lifetime tmax (decay with lifetime < tmax)
  //
  //fMaxLifetime = fConfig->GetDoubleDef("MaxLifetime", 1e-9);

  // Check whether the module is being run before or after the hadron
  // transport (intranuclear rescattering) module.
  //
  // If it is run before the hadron transport (and after the hadronization)
  // step it should decay only "unstable" particles (marked as hadrons in
  // the nucleus) which would typically decay within the time required to 
  // exit the nucleus - so, the algorithm wouldn't decay particles that 
  // have to be rescattered first. In case that the generated event is off
  // a free nucleon target, thi instance of the algorithm should do nothing.
  //  
  // If it is run after the hadon transport, then it should decay all the 
  // 'unstable' particles marked as 'present in the final state' and which 
  // should be decay before the event is passed to the detector particle
  // transport MC.
  //
  fRunBefHadroTransp = fConfig->GetBool("RunBeforeHadronTransport");

  // Allow user to specify a list of particles to be decayed
  //
  RgKeyList klist = gc->FindKeys("DecayParticleWithCode="); 
  RgKeyList::const_iterator kiter = klist.begin();
  for( ; kiter != klist.end(); ++kiter) { 
    RgKey key = *kiter;
    bool decay = gc->GetBool(key);
    vector<string> kv = utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdgc = atoi(kv[1].c_str());
    TParticlePDG * p = PDGLibrary::Instance()->Find(pdgc);
    if(decay) {
       LOG("ParticleDecayer", pDEBUG) 
            << "Configured to decay " <<  p->GetName();
       fParticlesToDecay.push_back(pdgc);
       vector <const DecayModelI *>::iterator diter = fDecayers->begin();
       for ( ; diter != fDecayers->end(); ++diter) {
           const DecayModelI * decayer = *diter;
           decayer->UnInhibitDecay(pdgc); 
       }// decayer
    }
    else {
       LOG("ParticleDecayer", pDEBUG) 
            << "Configured to inhibit decays for  " <<  p->GetName();
       fParticlesNotToDecay.push_back(pdgc);
       vector <const DecayModelI *>::iterator diter = fDecayers->begin();
       for ( ; diter != fDecayers->end(); ++diter) {
           const DecayModelI * decayer = *diter;
           decayer->InhibitDecay(pdgc); 
       }// decayer
    }// decay? 
  }// key iterator

  // Allow user to inhibit certain decay channels
  //
  klist = gc->FindKeys("InhibitDecay/"); 
  kiter = klist.begin();
  for( ; kiter != klist.end(); ++kiter) { 
    RgKey key = *kiter;
    if(gc->GetBool(key)) {
      string filtkey = utils::str::FilterString("InhibitDecay/", key);
      vector<string> kv = utils::str::Split(filtkey,",");
      assert(kv.size()==2);
      int pdgc = atoi(utils::str::FilterString("Particle=",kv[0]).c_str());
      int dc   = atoi(utils::str::FilterString("Channel=", kv[1]).c_str());
      TParticlePDG * p = PDGLibrary::Instance()->Find(pdgc);
      if(!p) continue;
      LOG("ParticleDecayer", pINFO) 
         << "Configured to inhibit " <<  p->GetName() 
         << "'s decay channel " << dc;
      vector <const DecayModelI *>::iterator diter = fDecayers->begin();
      for ( ; diter != fDecayers->end(); ++diter) {
         const DecayModelI * decayer = *diter;
         decayer->InhibitDecay(pdgc, p->DecayChannel(dc));
      }//decayer iterator
    }//val[key]=true?
  }//key iterator


  sort(fParticlesToDecay.begin(),    fParticlesToDecay.end());
  sort(fParticlesNotToDecay.begin(), fParticlesNotToDecay.end());

  // Print-out for only one of the two instances of this module
  if(!fRunBefHadroTransp) {
    LOG("ParticleDecayer", pNOTICE) 
       << "\nConfigured to decay: " << fParticlesToDecay
       << "\nConfigured to inhibit decays of: " << fParticlesNotToDecay
       << "\n";
  }
}
//___________________________________________________________________________
