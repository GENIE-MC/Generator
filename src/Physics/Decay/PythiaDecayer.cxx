//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <vector>
#include <cassert>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TDecayChannel.h>
#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif

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
#include "Physics/Decay/PythiaDecayer.h"

using std::vector;

using namespace genie;

// actual PYTHIA calls:
extern "C" void py1ent_(int *,  int *, double *, double *, double *);
extern "C" void pydecy_(int *);
//____________________________________________________________________________
PythiaDecayer::PythiaDecayer() :
Decayer("genie::PythiaDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::PythiaDecayer(string config) :
Decayer("genie::PythiaDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
PythiaDecayer::~PythiaDecayer()
{

}
//____________________________________________________________________________
void PythiaDecayer::ProcessEventRecord(GHepRecord * event) const
{
  LOG("ResonanceDecay", pINFO)
    << "Running PYTHIA6 particle decayer "
    << ((fRunBefHadroTransp) ? "*before*" : "*after*") << " FSI";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {
    ipos++;

    LOG("Pythia6Decay", pDEBUG) << "Checking: " << p->Name();

    int pdg_code = p->Pdg();
    GHepStatus_t status_code = p->Status();

    if(!this->IsHandled  (pdg_code)) continue;
    if(!this->ToBeDecayed(pdg_code, status_code)) continue;

    LOG("Pythia6Decay", pNOTICE)
          << "Decaying unstable particle: " << p->Name();

    bool decayed = this->Decay(ipos, event);
    assert(decayed); // handle this more graciously and throw an exception
  }

  LOG("Pythia6Decay", pNOTICE)
     << "Done finding & decaying unstable particles";
}
//____________________________________________________________________________
bool PythiaDecayer::Decay(int decay_particle_id, GHepRecord * event) const
{
  fWeight = 1.; // reset previous decay weight

  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return 0;

  // Get the particle 4-momentum, 4-position and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  TLorentzVector decay_particle_x4 = *(decay_particle->X4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Convert to PYTHIA6 particle code and check whether decay is inhibited
  int kc   = fPythia->Pycomp(decay_particle_pdg_code);
  int mdcy = fPythia->GetMDCY(kc, 1);
  if(mdcy == 0) {
    LOG("Pythia6Decay", pNOTICE)
       << (PDGLibrary::Instance())->Find(decay_particle_pdg_code)->GetName()
       << " decays are inhibited!";
    return false;
  }

  // Get sub of BRs and compute weight if decay channels were inhibited
  double sumbr = this->SumOfBranchingRatios(kc);
  if(sumbr <= 0) {
    LOG("Pythia6Decay", pWARN)
       << "The sum of enabled "
       << (PDGLibrary::Instance())->Find(decay_particle_pdg_code)->GetName()
       << " decay channel branching ratios is non-positive!";
    return false;
  }
  fWeight = 1./sumbr; // update weight to account for inhibited channels

  // Run PYTHIA6 decay
  int    ip    = 0;
  double E     = decay_particle_p4.Energy();
  double theta = decay_particle_p4.Theta();
  double phi   = decay_particle_p4.Phi();
  fPythia->SetMSTJ(22,1);
  py1ent_(&ip, &decay_particle_pdg_code, &E, &theta, &phi);

  // Get decay products
  fPythia->GetPrimaries();
  TClonesArray * impl = (TClonesArray *) fPythia->ImportParticles("All");
  if(!impl) {
    LOG("Pythia6Decay", pWARN)
      << "TPythia6::ImportParticles() returns a null list!";
    return false;
  }

  // Copy the PYTHIA6 container to the GENIE event record

  // Check whether the interaction is off a nuclear target or free nucleon
  // Depending on whether this module is run before or after the hadron
  // transport module it would affect the daughters status code
  GHepParticle * target_nucleus = event->TargetNucleus();
  bool in_nucleus = (target_nucleus!=0);

  TMCParticle * p = 0;
  TIter particle_iter(impl);
  while( (p = (TMCParticle *) particle_iter.Next()) ) {
    // Convert from TMCParticle to GHepParticle
    GHepParticle mcp = GHepParticle(
        p->GetKF(),                // pdg
        GHepStatus_t(p->GetKS()),  // status
        p->GetParent(),            // first parent
        0,                         // second parent
        p->GetFirstChild(),        // first daughter
        p->GetLastChild(),         // second daughter
        p->GetPx(),                // px
        p->GetPy(),                // py
        p->GetPz(),                // pz
        p->GetEnergy(),            // e
        p->GetVx(),                // x
        p->GetVy(),                // y
        p->GetVz(),                // z
        p->GetTime()               // t
    );

    if(mcp.Status()==kIStNucleonTarget) continue; // mother particle, already in GHEP

    int daughter_pdg_code = mcp.Pdg();
    SLOG("Pythia6Decay", pINFO)
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

  return true;
}
//____________________________________________________________________________
void PythiaDecayer::Initialize(void) const
{
  fPythia = TPythia6::Instance();
  fWeight = 1.;

  // sync GENIE/PYTHIA6 seeds
  RandomGen::Instance();
}
//____________________________________________________________________________
bool PythiaDecayer::IsHandled(int pdg_code) const
{
// does not handle requests to decay baryon resonances

 bool is_handled = (!utils::res::IsBaryonResonance(pdg_code));

 LOG("Pythia6Decay", pDEBUG)
    << "Can decay particle with PDG code = " << pdg_code
    << "? " << ((is_handled)? "Yes" : "No");

  return is_handled;
}
//____________________________________________________________________________
void PythiaDecayer::InhibitDecay(int pdg_code, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdg_code)) return;

  int kc = fPythia->Pycomp(pdg_code);

  if(!dc) {
    LOG("Pythia6Decay", pINFO)
       << "Switching OFF ALL decay channels for particle = " << pdg_code;
    fPythia->SetMDCY(kc, 1,0);
    return;
  }

  LOG("Pythia6Decay", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdg_code;

  int ichannel = this->FindPythiaDecayChannel(kc, dc);
  if(ichannel != -1) {
    fPythia->SetMDME(ichannel,1,0); // switch-off
  }
}
//____________________________________________________________________________
void PythiaDecayer::UnInhibitDecay(int pdg_code, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdg_code)) return;

  int kc = fPythia->Pycomp(pdg_code);

  if(!dc) {
    LOG("Pythia6Decay", pINFO)
      << "Switching ON all PYTHIA decay channels for particle = " << pdg_code;

    fPythia->SetMDCY(kc, 1,1);

    int first_channel = fPythia->GetMDCY(kc,2);
    int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

    for(int ichannel = first_channel;
            ichannel < last_channel; ichannel++) {
         fPythia->SetMDME(ichannel,1,1); // switch-on
    }
    return;
  }//!dc

  LOG("Pythia6Decay", pINFO)
     << "Switching OFF decay channel = " << dc->Number()
     << " for particle = " << pdg_code;

  int ichannel = this->FindPythiaDecayChannel(kc, dc);
  if(ichannel != -1) {
    fPythia->SetMDME(ichannel,1,1); // switch-on
  }
}
//____________________________________________________________________________
double PythiaDecayer::SumOfBranchingRatios(int kc) const
{
// Sum of branching ratios for enabled channels
//
  double sumbr=0.;

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  bool has_inhibited_channels=false;

  // loop over pythia decay channels
  for(int ichannel = first_channel;
          ichannel < last_channel; ichannel++) {

     bool enabled = (fPythia->GetMDME(ichannel,1) == 1);
     if (!enabled) {
       has_inhibited_channels = true;
     } else {
       sumbr += fPythia->GetBRAT(ichannel);
     }
  }

  if(!has_inhibited_channels) return 1.;
  LOG("Pythia6Decay", pINFO) << "Sum{BR} = " << sumbr;

  return sumbr;
}
//____________________________________________________________________________
int PythiaDecayer::FindPythiaDecayChannel(int kc, TDecayChannel* dc) const
{
  if(!dc) return -1;

  int first_channel = fPythia->GetMDCY(kc,2);
  int last_channel  = fPythia->GetMDCY(kc,2) + fPythia->GetMDCY(kc,3) - 1;

  bool found_match = false;

  // loop over pythia decay channels
  for(int ichannel = first_channel;
          ichannel < last_channel; ichannel++) {

     // does the  current pythia channel matches the input TDecayChannel?
     LOG("Pythia6Decay", pINFO)
         << "\nComparing PYTHIA's channel = " << ichannel
         << " with TDecayChannel = " << dc->Number();

     found_match = this->MatchDecayChannels(ichannel,dc);
     if(found_match) {
         LOG("Pythia6Decay", pNOTICE)
            << " ** TDecayChannel id = " << dc->Number()
            << " corresponds to PYTHIA6 channel id = " << ichannel;
         return ichannel;
     }//match?
  }//loop pythia decay ch.

  LOG("Pythia6Decay", pWARN)
     << " ** No PYTHIA6 decay channel match found for "
     << "TDecayChannel id = " << dc->Number();

  return -1;
}
//____________________________________________________________________________
bool PythiaDecayer::MatchDecayChannels(int ichannel, TDecayChannel* dc) const
{
  // num. of daughters in the input TDecayChannel & the input PYTHIA ichannel
  int nd = dc->NDaughters();

  int py_nd = 0;
  for (int i = 1; i <= 5; i++) {
     if(fPythia->GetKFDP(ichannel,i) != 0) py_nd++;
  }

  LOG("Pythia6Decay", pDEBUG)
    << "NDaughters: PYTHIA = " << py_nd << ", ROOT's TDecayChannel = " << nd;

  if(nd != py_nd) return false;

  //
  // if the two channels have the same num. of daughters, then compare them
  //

  // store decay daughters for the input TDecayChannel
  vector<int> dc_daughter(nd);
  int k=0;
  for( ; k < nd; k++) {
     dc_daughter[k] = dc->DaughterPdgCode(k);
  }
  // store decay daughters for the input PYTHIA's ichannel
  vector<int> py_daughter(nd);
  k=0;
  for(int i = 1; i <= 5; i++) {
     if(fPythia->GetKFDP(ichannel,i) == 0) continue;
     py_daughter[k] = fPythia->GetKFDP(ichannel,i);
     k++;
  }

  // sort both daughter lists
  sort( dc_daughter.begin(), dc_daughter.end() );
  sort( py_daughter.begin(), py_daughter.end() );

  // compare
  for(int i = 0; i < nd; i++) {
    if(dc_daughter[i] != py_daughter[i]) return false;
  }
  return true;
}
//____________________________________________________________________________
