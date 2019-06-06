//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <algorithm>
#include <sstream>

#include <TParticlePDG.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"
#include "Physics/Decay/Decayer.h"

using std::count;
using std::ostringstream;

using namespace genie;
//___________________________________________________________________________
Decayer::Decayer() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
Decayer::Decayer(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
Decayer::Decayer(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
Decayer::~Decayer()
{

}
//___________________________________________________________________________
bool Decayer::ToBeDecayed(int pdg_code, GHepStatus_t status_code) const
{
  // Check whether it is "unstable" (definition can vary)

  bool is_unstable = this->IsUnstable(pdg_code);

  LOG("Decay", pDEBUG)
    << "Particle is unstable? "
    << ((is_unstable) ? "Yes" : "No");

  if(!is_unstable) return false;

  // Check whether the given unstable particle
  // has the appropriate status code to be decayed

  bool to_be_decayed = false;

  if(fRunBefHadroTransp) {
    to_be_decayed =
      (status_code == kIStHadronInTheNucleus    ||
       status_code == kIStPreDecayResonantState ||
       status_code == kIStStableFinalState);
  }
  else {
    to_be_decayed =
      (status_code == kIStStableFinalState);
  }

  LOG("Decay", pDEBUG)
      << "Particle to be decayed "
      << "[" << ((fRunBefHadroTransp) ? "Before" : "After") << " FSI]? "
      << ((to_be_decayed) ? "Yes" : "No");

  return to_be_decayed;
}
//___________________________________________________________________________
bool Decayer::IsUnstable(int pdg_code) const
{
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
void Decayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void Decayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void Decayer::LoadConfig(void)
{
  // Get the specified maximum lifetime tmax (decay with lifetime < tmax)
  //
  //fMaxLifetime = fConfig->GetDoubleDef("MaxLifetime", 1e-9);

  // Check whether to generate weighted or unweighted particle decays
  fGenerateWeighted = false ;
  //this->GetParam("GenerateWeighted", fGenerateWeighted, false);

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
  this->GetParam("RunBeforeHadronTransport", fRunBefHadroTransp) ;

  // Allow user to specify a list of particles to be decayed
  //
  RgKeyList klist = GetConfig().FindKeys("DecayParticleWithCode=");
  RgKeyList::const_iterator kiter = klist.begin();
  for( ; kiter != klist.end(); ++kiter) {
    RgKey key = *kiter;
    bool decay = GetConfig().GetBool(key);
    vector<string> kv = utils::str::Split(key,"=");
    assert(kv.size()==2);
    int pdgc = atoi(kv[1].c_str());
    TParticlePDG * p = PDGLibrary::Instance()->Find(pdgc);
    if(decay) {
       LOG("Decay", pDEBUG)
            << "Configured to decay " <<  p->GetName();
       fParticlesToDecay.push_back(pdgc);
       this->UnInhibitDecay(pdgc);
    }
    else {
       LOG("Decay", pDEBUG)
            << "Configured to inhibit decays for  " <<  p->GetName();
       fParticlesNotToDecay.push_back(pdgc);
       this->InhibitDecay(pdgc);
    }// decay?
  }// key iterator

  // Allow user to inhibit certain decay channels
  //
  klist = GetConfig().FindKeys("InhibitDecay/");
  kiter = klist.begin();
  for( ; kiter != klist.end(); ++kiter) {
    RgKey key = *kiter;
    if(GetConfig().GetBool(key)) {
      string filtkey = utils::str::FilterString("InhibitDecay/", key);
      vector<string> kv = utils::str::Split(filtkey,",");
      assert(kv.size()==2);
      int pdgc = atoi(utils::str::FilterString("Particle=",kv[0]).c_str());
      int dc   = atoi(utils::str::FilterString("Channel=", kv[1]).c_str());
      TParticlePDG * p = PDGLibrary::Instance()->Find(pdgc);
      if(!p) continue;
      LOG("Decay", pINFO)
         << "Configured to inhibit " <<  p->GetName()
         << "'s decay channel " << dc;
      this->InhibitDecay(pdgc, p->DecayChannel(dc));
    }//val[key]=true?
  }//key iterator


  sort(fParticlesToDecay.begin(),    fParticlesToDecay.end());
  sort(fParticlesNotToDecay.begin(), fParticlesNotToDecay.end());

  // Print-out for only one of the two instances of this module
  if(!fRunBefHadroTransp) {
    LOG("Decay", pNOTICE)
       << "\nConfigured to decay: " << fParticlesToDecay
       << "\nConfigured to inhibit decays of: " << fParticlesNotToDecay
       << "\n";
  }
}
//___________________________________________________________________________
