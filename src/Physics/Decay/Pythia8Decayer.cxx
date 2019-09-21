//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
         Queen Mary University of London
*/
//____________________________________________________________________________

#include <vector>
#include <cassert>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TDecayChannel.h>
#include <RVersion.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/Decay/Pythia8Decayer.h"

using std::vector;

using namespace genie;

//____________________________________________________________________________
Pythia8Decayer::Pythia8Decayer() :
Decayer("genie::Pythia8Decayer")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Decayer::Pythia8Decayer(string config) :
Decayer("genie::Pythia8Decayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Decayer::~Pythia8Decayer()
{

}
//____________________________________________________________________________
void Pythia8Decayer::ProcessEventRecord(GHepRecord * event) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  LOG("ResonanceDecay", pINFO)
    << "Running PYTHIA8 particle decayer "
    << ((fRunBefHadroTransp) ? "*before*" : "*after*") << " FSI";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {
    ipos++;

    LOG("Pythia8Decay", pDEBUG) << "Checking: " << p->Name();

    int pdg_code = p->Pdg();
    GHepStatus_t status_code = p->Status();

    if(!this->IsHandled  (pdg_code)) continue;
    if(!this->ToBeDecayed(pdg_code, status_code)) continue;

    LOG("Pythia8Decay", pNOTICE)
          << "Decaying unstable particle: " << p->Name();

    bool decayed = this->Decay(ipos, event);
    assert(decayed); // handle this more graciously and throw an exception
  }

  LOG("Pythia8Decay", pNOTICE)
     << "Done finding & decaying unstable particles";
#endif
}
//____________________________________________________________________________
bool Pythia8Decayer::Decay(int decay_particle_id, GHepRecord * event) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fWeight = 1.; // reset previous decay weight

  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return 0;

  // Get the particle 4-momentum, 4-position and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  TLorentzVector decay_particle_x4 = *(decay_particle->X4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Convert to PYTHIA8 particle code and check whether decay is inhibited
  bool md = fPythia->Pythia8()->particleData.canDecay(decay_particle_pdg_code);
  if(not md) {
    LOG("Pythia8Decay", pNOTICE)
       << (PDGLibrary::Instance())->Find(decay_particle_pdg_code)->GetName()
       << " decays are inhibited!";
    return false;
  }

  // Get sub of BRs and compute weight if decay channels were inhibited
  double sumbr = this->SumOfBranchingRatios(decay_particle_pdg_code);
  if(sumbr <= 0) {
    LOG("Pythia8Decay", pWARN)
       << "The sum of enabled "
       << (PDGLibrary::Instance())->Find(decay_particle_pdg_code)->GetName()
       << " decay channel branching ratios is non-positive!";
    return false;
  }
  fWeight = 1./sumbr; // update weight to account for inhibited channels

  // Run PYTHIA8 decay
  double E     = decay_particle_p4.Energy();
  double M     = fPythia->Pythia8()->particleData.m0(decay_particle_pdg_code);
  double pz    = Pythia8::sqrtpos(E*E - M*M);

  fPythia->Pythia8()->event.reset();

  fPythia->Pythia8()->event.append(
      decay_particle_pdg_code, 11, 0, 0, 0., 0., pz, E, M
  );
  fPythia->Pythia8()->next();

  // List event information
  fPythia->Pythia8()->event.list();
  fPythia->Pythia8()->stat();

  //-- get decay products
  Pythia8::Event &fEvent = fPythia->Pythia8()->event;
  int numpart = fEvent.size();

  // Check whether the interaction is off a nuclear target or free nucleon
  // Depending on whether this module is run before or after the hadron
  // transport module it would affect the daughters status code
  GHepParticle * target_nucleus = event->TargetNucleus();
  bool in_nucleus = (target_nucleus!=0);

  int ioff = 0;
  if (fEvent[0].id() == 90) ioff = -1;

  for (int i = 1; i < numpart; ++i) {
    /*
     * Convert Pythia8 status code to Pythia6
     * Decayed/fragmented particle has a pytahi6 code of 11 (kIStNucleonTarget)
     * Final state particles have a negative pythia8 code and a pythia6 code of 1 (kIStStableFinalState)
     */
    GHepStatus_t gStatus = (fEvent[i].status()>0) ? kIStStableFinalState : kIStNucleonTarget;
    GHepParticle mcp = GHepParticle(
            fEvent[i].id(),                // pdg
            gStatus,                       // status
            fEvent[i].mother1()   + ioff,  // first parent
            fEvent[i].mother2()   + ioff,  // second parent
            fEvent[i].daughter1() + ioff,  // first daughter
            fEvent[i].daughter2() + ioff,  // second daughter
            fEvent[i].px(),                // px [GeV/c]
            fEvent[i].py(),                // py [GeV/c]
            fEvent[i].pz(),                // pz [GeV/c]
            fEvent[i].e(),                 // e  [GeV]
            fEvent[i].xProd(),             // x  [mm]
            fEvent[i].yProd(),             // y  [mm]
            fEvent[i].zProd(),             // z  [mm]
            fEvent[i].tProd()              // t  [mm/c]
    );

    if(mcp.Status()==kIStNucleonTarget) continue; // mother particle, already in GHEP

    int daughter_pdg_code = mcp.Pdg();
    SLOG("Pythia8Decay", pINFO)
       << "Adding daughter particle wit PDG code = "
       << daughter_pdg_code << ", m = " << mcp.Mass()
       << " GeV, E = " << mcp.Energy() << " GeV)";

    bool is_hadron = pdg::IsHadron(daughter_pdg_code);
    bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

    GHepStatus_t daughter_status_code = (hadron_in_nuc) ?
         kIStHadronInTheNucleus : kIStStableFinalState;

    TLorentzVector daughter_p4(
       mcp.Px(),mcp.Py(),mcp.Pz(),mcp.Energy());
    event->AddParticle(
       daughter_pdg_code, daughter_status_code,
       decay_particle_id,-1,-1,-1,
       daughter_p4, decay_particle_x4);
  }

  // Update the event weight for each weighted particle decay
  double weight = event->Weight() * fWeight;
  event->SetWeight(weight);

  // Mark input particle as a 'decayed state' & add its daughter links
  decay_particle->SetStatus(kIStDecayedState);

#endif
  return true;
}
//____________________________________________________________________________
void Pythia8Decayer::Initialize(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = Pythia8Singleton::Instance();
  fWeight = 1.;
  fPythia->Pythia8()->readString("ProcessLevel:all = off");
  fPythia->Pythia8()->readString("Print:quiet      = on");

  // sync GENIE/PYTHIA8 seeds
  RandomGen::Instance();
  fPythia->Pythia8()->init();
#endif
}
//____________________________________________________________________________
bool Pythia8Decayer::IsHandled(int pdg_code) const
{
// does not handle requests to decay baryon resonances

 bool is_handled = (!utils::res::IsBaryonResonance(pdg_code));

 LOG("Pythia8Decay", pDEBUG)
    << "Can decay particle with PDG code = " << pdg_code
    << "? " << ((is_handled)? "Yes" : "No");

  return is_handled;
}
//____________________________________________________________________________
void Pythia8Decayer::InhibitDecay(int pdg_code, TDecayChannel * dc) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  if(! this->IsHandled(pdg_code)) return;

  if(!dc) {
    LOG("Pythia8Decay", pINFO)
       << "Switching OFF ALL decay channels for particle = " << pdg_code;
    fPythia->Pythia8()->particleData.mayDecay(pdg_code, false);
    return;
  }

  LOG("Pythia8Decay", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdg_code;

  int ichannel = this->FindPythiaDecayChannel(pdg_code, dc);
  if(ichannel != -1) {
    Pythia8::ParticleDataEntry * fPDE =
      fPythia->Pythia8()->particleData.particleDataEntryPtr(pdg_code);
    fPDE->channel(ichannel).onMode(0); // switch-off
  }
#endif
}
//____________________________________________________________________________
void Pythia8Decayer::UnInhibitDecay(int pdg_code, TDecayChannel * dc) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  if(! this->IsHandled(pdg_code)) return;

  Pythia8::ParticleDataEntry * fPDE =
    fPythia->Pythia8()->particleData.particleDataEntryPtr(pdg_code);

  if(!dc) {
    LOG("Pythia8Decay", pINFO)
      << "Switching ON all PYTHIA decay channels for particle = " << pdg_code;

    fPythia->Pythia8()->particleData.mayDecay(pdg_code, true);

    // loop over pythia decay channels
    int size_channels = fPDE->sizeChannels();
    for (int ichannel = 0; ichannel < size_channels; ++ichannel) {
      fPDE->channel(ichannel).onMode(1); // switch-on
    }

    return;
  }//!dc

  LOG("Pythia8Decay", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdg_code;

  int ichannel = this->FindPythiaDecayChannel(pdg_code, dc);
  if(ichannel != -1) {
    fPDE->channel(ichannel).onMode(1); // switch-on
  }
#endif
}
//____________________________________________________________________________
double Pythia8Decayer::SumOfBranchingRatios(int pdg_code) const
{
// Sum of branching ratios for enabled channels
//
  double sumbr=0.;
#ifdef __GENIE_PYTHIA8_ENABLED__

  bool has_inhibited_channels=false;

  Pythia8::ParticleDataEntry * fPDE =
    fPythia->Pythia8()->particleData.particleDataEntryPtr(pdg_code);
  int size_channels = fPDE->sizeChannels();

  // loop over pythia decay channels
  for (int ichannel = 0; ichannel < size_channels; ++ichannel) {
     bool enabled = (fPDE->channel(ichannel).onMode() == 1);
     if (!enabled) {
       has_inhibited_channels = true;
     } else {
       sumbr += fPDE->channel(ichannel).bRatio();
     }
  }

  if(!has_inhibited_channels) return 1.;
  LOG("Pythia8Decay", pINFO) << "Sum{BR} = " << sumbr;

#endif
  return sumbr;
}
//____________________________________________________________________________
int Pythia8Decayer::FindPythiaDecayChannel(int pdg_code, TDecayChannel* dc) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  if(!dc) return -1;

  bool found_match = false;

  Pythia8::ParticleDataEntry * fPDE =
    fPythia->Pythia8()->particleData.particleDataEntryPtr(pdg_code);
  int size_channels = fPDE->sizeChannels();

  // loop over pythia decay channels
  for (int ichannel = 0; ichannel < size_channels; ++ichannel) {

     // does the  current pythia channel match the input TDecayChannel?
     LOG("Pythia8Decay", pINFO)
         << "\nComparing PYTHIA's channel = " << ichannel
         << " with TDecayChannel = " << dc->Number();

     found_match = this->MatchDecayChannels(pdg_code, ichannel, dc);
     if(found_match) {
         LOG("Pythia8Decay", pNOTICE)
            << " ** TDecayChannel id = " << dc->Number()
            << " corresponds to PYTHIA8 channel id = " << ichannel;
         return ichannel;
     }//match?
  }//loop pythia decay ch.

  LOG("Pythia8Decay", pWARN)
     << " ** No PYTHIA8 decay channel match found for "
     << "TDecayChannel id = " << dc->Number();

#endif
  return -1;
}
//____________________________________________________________________________
bool Pythia8Decayer::MatchDecayChannels(
    int pdg_code, int ichannel, TDecayChannel* dc) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  // num. of daughters in the input TDecayChannel & the input PYTHIA ichannel
  int nd = dc->NDaughters();

  Pythia8::ParticleDataEntry * fPDE =
    fPythia->Pythia8()->particleData.particleDataEntryPtr(pdg_code);

  int py_nd = 0;
  for (int i = 1; i <= 8; i++) {
    if(fPDE->channel(ichannel).product(i) != 0) py_nd++;
  }

  LOG("Pythia8Decay", pDEBUG)
    << "NDaughters: PYTHIA = " << py_nd << ", ROOT's TDecayChannel = " << nd;

  if(nd != py_nd) return false;

  //
  // if the two channels have the same num. of daughters, then compare them
  //

  // store decay daughters for the input TDecayChannel
  vector<int> dc_daughter(nd);
  for (int i = 0; i < nd; i++) {
     dc_daughter[i] = dc->DaughterPdgCode(i);
  }
  // store decay daughters for the input PYTHIA's ichannel
  vector<int> py_daughter(nd);
  for (int i = 0; i < py_nd; i++) {
      py_daughter[i] = fPDE->channel(ichannel).product(i);
  }

  // sort both daughter lists
  sort( dc_daughter.begin(), dc_daughter.end() );
  sort( py_daughter.begin(), py_daughter.end() );

  // compare
  for(int i = 0; i < nd; i++) {
    if(dc_daughter[i] != py_daughter[i]) return false;
  }
#endif
  return true;
}
//____________________________________________________________________________
