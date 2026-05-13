//____________________________________________________________________________
/*
 Copyright (c) 2003-2026, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <vector>
#include <cassert>
#include <algorithm>

#include <iostream>
#include <sstream>
#include <iomanip>

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
#include "Physics/Decay/Pythia8Decayer2023.h"

using std::vector;

using namespace genie;

//____________________________________________________________________________
Pythia8Decayer2023::Pythia8Decayer2023() :
Decayer("genie::Pythia8Decayer2023")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Decayer2023::Pythia8Decayer2023(string config) :
Decayer("genie::Pythia8Decayer2023", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Decayer2023::~Pythia8Decayer2023()
{

}
//____________________________________________________________________________
void Pythia8Decayer2023::ProcessEventRecord(GHepRecord * event) const
{
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
}
//____________________________________________________________________________
bool Pythia8Decayer2023::Decay(int decay_particle_id, GHepRecord * event) const
{
  fWeight = 1.; // reset previous decay weight

#ifdef __GENIE_PYTHIA8_ENABLED__

  // Get the one-and-only instance of Pythia8 that GENIE will use
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return 0;

  // Check whether the interaction is off a nuclear target or free nucleon
  // Depending on whether this module is run before or after the hadron
  // transport module it would affect the daughters status code
  GHepParticle * target_nucleus = event->TargetNucleus();
  bool in_nucleus = (target_nucleus!=0);

  // Get the particle 4-momentum, 4-position and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  TLorentzVector decay_particle_x4 = *(decay_particle->X4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  TVector3 polz;
  decay_particle->GetPolarization(polz);

  gPythia->event.reset();

  // check if pdgid is consistent with Pythia8 particle table
  //
  if ( !gPythia->particleData.findParticle(decay_particle_pdg_code) ) {
    LOG("Pythia8Decay", pWARN)
      << " can NOT find pdgid = " << decay_particle_pdg_code
      << " in Pythia8::ParticleData";
      return false;
    }

  if ( !gPythia->particleData.canDecay(decay_particle_pdg_code) ) {
    LOG("Pythia8Decay", pWARN)
      << " Particle of pdgid = " << decay_particle_pdg_code
      << " can NOT be decayed by Pythia8";
      return false;
  }

  // NOTE: Energy should be in GeV

  gPythia->event.append( decay_particle_pdg_code, 1, 0, 0,
                         decay_particle_p4.Px(),
                         decay_particle_p4.Py(),
                         decay_particle_p4.Pz(),
                         decay_particle_p4.E(),
                         decay_particle_p4.M() );

  // specify polarization, if any

  // NOTE: while in Py8 polarization is a double variable ,
  //       in reality it's expected to be -1, 0., or 1 in case of "external" tau's,
  //       similar to LHA SPINUP; see Particle Decays, Hadron and Tau Decays in docs at
  //       https://pythia.org/manuals/pythia8305/Welcome.html
  //       so it's not able to handle anything like 0.99, thus we're rounding off
  //  gPythia->event.back().pol( round( std::cos( track.GetPolarization().angle( track.GetMomentumDirection() ) ) ) );
  double a = decay_particle_p4.Vect().Angle(polz);
  gPythia->event.back().pol( round( std::cos(a) ) );

  int npart_before_decay = gPythia->event.size();

  gPythia->next();

  //gPythia->event.list();
  //gPythia->stat();

  int npart_after_decay = gPythia->event.size();

  LOG("Pythia8Decay", pINFO) << "before " << npart_before_decay << " after " << npart_after_decay;

  // the values of space coordinates from pythia are in mm.
  // our conventions want it in fm
  constexpr double space_scale = units::mm / units::fm ;

  //  the values of time coordinate from pythia is in mm/c.
  // our conventions want it in ys
  constexpr double time_scale = 1e21 * units::m / units::s ;


  // create & fill up decay products
  //
  for ( int ip=npart_before_decay; ip<npart_after_decay; ++ip )
    {

      // only select final state decay products (direct or via subsequent decays);
      // skip all others
      //
      // NOTE: in general, final state decays products will have
      //       positive status code between 91 and 99
      //       (in case such information could be of interest in the future)
      //
      if ( gPythia->event[ip].status() < 0 ) continue;

      Pythia8::Event &fEvent = gPythia->event;

      int daughter_pdg_code   = fEvent[ip].id();

      int mother1 = decay_particle_id;
      int mother2 = -1;

      int daughter1 = -1;
      int daughter2 = -1;

      bool is_hadron = pdg::IsHadron(daughter_pdg_code);
      bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

      GHepStatus_t daughter_status_code = (hadron_in_nuc) ?
        kIStHadronInTheNucleus : kIStStableFinalState;

      GHepParticle mcp = GHepParticle(
                                      daughter_pdg_code, // pdg
                                      daughter_status_code, // status
                                      mother1,           // first parent
                                      mother2,           // second parent
                                      daughter1,         // first daughter
                                      daughter2,         // second daughter
                                      fEvent[ip].px(),   // px
                                      fEvent[ip].py(),   // py
                                      fEvent[ip].pz(),   // pz
                                      fEvent[ip].e(),    // e
                                      fEvent[ip].xProd() * space_scale + decay_particle_x4.X(), // x
                                      fEvent[ip].yProd() * space_scale + decay_particle_x4.Y(), // y
                                      fEvent[ip].zProd() * space_scale + decay_particle_x4.Z(), // z
                                      fEvent[ip].tProd() * time_scale  + decay_particle_x4.T()  // t
                                      );

    SLOG("Pythia8Decay", pINFO)
       << "Adding daughter particle with PDG code = "
       << daughter_pdg_code << ", m = " << mcp.Mass()
       << " GeV, E = " << mcp.Energy() << " GeV)";

    event->AddParticle(mcp);

    }

  // Update the event weight for each weighted particle decay
  double weight = event->Weight() * fWeight;
  event->SetWeight(weight);

  // Mark input particle as a 'decayed state' & add its daughter links
  decay_particle->SetStatus(kIStDecayedState);

  return true;
#else
  (void) decay_particle_id;
  (void) event;
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);

  return false;
#endif

}
//____________________________________________________________________________
void Pythia8Decayer2023::Initialize(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fWeight = 1.;

  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();
  gPythia->readString("ProcessLevel:all = off");
  gPythia->readString("ProcessLevel:resonanceDecays=on");

  gPythia->readString("Print:quiet = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  gPythia->readString("Random::setSeed = on");
  gPythia->settings.mode("Random:seed",seed);
  LOG("Pythia8Decay", pINFO)
    << "PYTHIA8 seed = " << gPythia->settings.mode("Random:seed");

  gPythia->init();

  // shut off decays of pi0's as we want Geant4 to handle them
  // if other immediate decay products should be handled by Geant4,
  // their respective decay modes should be shut off as well
  //
  gPythia->readString("111:onMode = off");

#else
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
bool Pythia8Decayer2023::IsHandled(int pdg_code) const
{
// does not handle requests to decay baryon resonances

 bool is_handled = (!utils::res::IsBaryonResonance(pdg_code));

 LOG("Pythia8Decay", pDEBUG)
    << "Can decay particle with PDG code = " << pdg_code
    << "? " << ((is_handled)? "Yes" : "No");

  return is_handled;
}
//____________________________________________________________________________
void Pythia8Decayer2023::InhibitDecay(int pdg_code, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdg_code)) return;

  LOG("Pythia8Decay",pINFO)
    << "InhibitDecay for pdgc " << pdg_code << " PDGLibrary channel "
    << (dc?dc->Number():-999) << " [ -999 = all ]";

#ifdef __GENIE_PYTHIA8_ENABLED__

  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();
  bool known = gPythia->particleData.isParticle(pdg_code);
  if ( ! known ) {
    LOG("Pythia8Decay", pERROR)
      << "Can not switch off decays of pdgc " << pdg_code
      << " as it is UNKNOWN to Pythia8";
    return;
  }

  auto pdentry = gPythia->particleData.particleDataEntryPtr(pdg_code);

  int ifirst_chan = 0;
  int ilast_chan  = pdentry->sizeChannels() - 1;
  if (ilast_chan < 0) {
    LOG("Pythia8Decay", pDEBUG)
      << "No available decay channels for pdgc " << pdg_code;
    return;
  }

  if (dc) {
    // find the matching channel set ifirst_chan,ilast_chan to that
    LOG("Pythia8Decay",pDEBUG)
      << "Need to map TDatabasePDG channel to Pythia8 channel for " << pdg_code;
    int ichannel = FindPythiaDecayChannel(pdg_code,dc);
    if (ichannel < 0) {
      LOG("Pythia8Decay",pWARN)
        << "Could not map TDecayChannel " << dc->Number() << " to a Pythia8 channel";
      return;
    }
    LOG("Pythia8Decay",pINFO)
      << "TDatabasePDG channel " << dc->Number() << " mapped to Pythia8 " << ichannel;
    // toggle just the one channel
    ifirst_chan = ichannel;
    ilast_chan  = ichannel;
  }

  LOG("Pythia8Decay",pINFO)
    << "Inhibiting Pythia8 channels [" << ifirst_chan << "," << ilast_chan << "] of pdgc " << pdg_code;
  //std::cout << "Before inhibiting channels " << ifirst_chan << "," << ilast_chan
  //          << " of " << pdg_code << std::endl;
  //gPythia->particleData.list(pdg_code);

  for (int ichan=ifirst_chan; ichan<=ilast_chan; ++ichan) {
    int onMode = pdentry->channel(ichan).onMode();
    int onMode_original = onMode;
    bool is_particle = (pdg_code>0);
    // 4 possible modes:  0=off particle & antiparticle; 1=on both; 2=on particle only; 3=on antiparticle only
    // we're trying to turn _off_ the mode
    switch ( onMode ) {
    case 0: // off for particles & antiparticles
      break; // already off
    case 1: // on for both particle & antiparticle
      onMode = (is_particle ? 3 : 2);
      break;
    case 2: // on for particle only
      onMode = (is_particle ? 0 : 2);
      break;
    case 3: // on for antiparticle only
      onMode = (is_particle ? 3 : 0);
      break;
    }
    LOG("Pythia8Decay", pINFO) << " change py8 ichan " << ichan << ": " << onMode_original << " to " << onMode;
    pdentry->channel(ichan).onMode(onMode);
  }

  //std::cout << "*** After inhibiting channels " << ifirst_chan << "," << ilast_chan
  //          << " of " << pdg_code << std::endl;
  //gPythia->particleData.list(pdg_code);

#else
  (void) dc;
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif

}
//____________________________________________________________________________
void Pythia8Decayer2023::UnInhibitDecay(int pdg_code, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdg_code)) return;

  LOG("Pythia8Decay",pINFO)
    << "UnInhibitDecay for pdgc " << pdg_code << " PDGLibrary channel "
    << (dc?dc->Number():-999) << " [ -999 = all ]";

#ifdef __GENIE_PYTHIA8_ENABLED__

  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();
  bool known = gPythia->particleData.isParticle(pdg_code);
  if ( ! known ) {
    LOG("Pythia8Decay", pERROR)
      << "Can not switch on decays of pdgc " << pdg_code
      << " as it is UNKNOWN to Pythia8";
    return;
  }

  auto pdentry = gPythia->particleData.particleDataEntryPtr(pdg_code);

  int ifirst_chan = 0;
  int ilast_chan  = pdentry->sizeChannels() - 1;
  if (ilast_chan < 0) {
    LOG("Pythia8Decay", pDEBUG)
      << "No available decay channels for pdgc " << pdg_code;
    return;
  }

  if (dc) {
    // find the matching channel set ifirst_chan,ilast_chan to that
    LOG("Pythia8Decay",pDEBUG)
      << "Need to map TDatabasePDG channel to Pythia8 channel for " << pdg_code << " (" << dc->Number() << ")";
    int ichannel = FindPythiaDecayChannel(pdg_code,dc);
    if (ichannel < 0) {
      LOG("Pythia8Decay",pWARN)
        << "Could not map TDecayChannel " << dc->Number() << " to a Pythia8 channel";
      return;
    }
    LOG("Pythia8Decay",pINFO)
      << "TDatabasePDG channel " << dc->Number() << " mapped to Pythia8 " << ichannel;
    // toggle just the one channel
    ifirst_chan = ichannel;
    ilast_chan  = ichannel;
  }

  LOG("Pythia8Decay",pINFO)
    << "Unihibiting Pythia8 channels [" << ifirst_chan << "," << ilast_chan << "] of pdgc " << pdg_code;
  //std::cout << "Before uninhibiting channels " << ifirst_chan << "," << ilast_chan
  //          << " of " << pdg_code << std::endl;
  //gPythia->particleData.list(pdg_code);

  for (int ichan=ifirst_chan; ichan<=ilast_chan; ++ichan) {
    int onMode = pdentry->channel(ichan).onMode();
    int onMode_original = onMode;
    bool is_particle = (pdg_code>0);
    // 4 possible modes:  0=off particle & antiparticle; 1=on both; 2=on particle only; 3=on antiparticle only
    // we're trying to turn _on_ the mode
    switch ( onMode ) {
    case 0: // off for particles & antiparticles
      onMode = (is_particle ? 2 : 3);
      break;
    case 1: // on for both particle & antiparticle
      break; // already on
    case 2: // on for particle only
      onMode = (is_particle ? 2 : 1);
      break;
    case 3: // on for antiparticle only
      onMode = (is_particle ? 1 : 3);
      break;
    }
    LOG("Pythia8Decay", pINFO) << " change py8 ichan " << ichan << ": " << onMode_original << " to " << onMode;
    pdentry->channel(ichan).onMode(onMode);
  }

  //std::cout << "After uninhibiting channels " << ifirst_chan << "," << ilast_chan
  //          << " of " << pdg_code << std::endl;
  //gPythia->particleData.list(pdg_code);

#else
  (void) dc;
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
double Pythia8Decayer2023::SumOfBranchingRatios(int /* kc */ ) const
{
// Sum of branching ratios for enabled channels
//
  double sumbr=0.;

#ifdef __GENIE_PYTHIA8_ENABLED__

  LOG("Pythia8Decay", pFATAL)
    << "Pythia8Decayer2023::SumOfBranchingRatios() not yet implemented";
  gAbortingInErr = true;
  std::exit(1);

  LOG("Pythia8Decay", pINFO) << "Sum{BR} = " << sumbr;

#else
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif

  return sumbr;
}
//____________________________________________________________________________
int Pythia8Decayer2023::FindPythiaDecayChannel(int pdgc, TDecayChannel* dc) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  LOG("Pythia8Decay",pDEBUG)
    << "looking for " << pdgc << " channel " << dc->Number();

  // Get the one-and-only instance of Pythia8 that GENIE will use
  Pythia8::Pythia* gPythia = genie::Pythia8Singleton::Instance()->Pythia8();
  if ( ! gPythia->particleData.isParticle(pdgc) ) {
    LOG("Pythia8Decay",pWARN) << "Pythia8 doesn't know about pdg " << pdgc;
    return -1;
  }

  // gather info from TDatabasePDG channel; get sorted list of daughters to match
  int nd = dc->NDaughters();
  vector<int> dc_daughter(nd);
  for (int i=0; i < nd; ++i) dc_daughter[i] = dc->DaughterPdgCode(i);
  sort( dc_daughter.begin(), dc_daughter.end() );

  auto pdentry = gPythia->particleData.particleDataEntryPtr(pdgc);
  // loop over Pythia8 channels looking for a match
  int nChan = pdentry->sizeChannels();

  LOG("Pythia8Decay",pDEBUG) << " dc has nd " << nd << " daughters "
                             << " pythia8 has nChan " << nChan;

  for (int iChan=0; iChan < nChan; ++iChan) {
    auto p8dkchan = pdentry->channel(iChan);

    int nmult = p8dkchan.multiplicity();
    // must have same number of daughters
    if (nd != nmult) {
      LOG("Pythia8Decay",pDEBUG) << "iChan " << iChan << " " << nd << " != " << nmult;
      continue;
    }

    vector<int> py_daughter(nmult);
    for (int j=0; j < nmult; ++j) {
      // pythia8 decay list is for particles, not antiparticles
      // if we're getting the daughters of an antiparticle, then reverse them
      int py_dau = p8dkchan.product(j);
      if (pdgc<0) py_dau = gPythia->particleData.antiId(py_dau);
      py_daughter[j] = py_dau;
    }
    sort ( py_daughter.begin(), py_daughter.end() );

    // now compare the lists
    bool match = true;
    for (int k=0; k < nd; ++k) {
      LOG("Pythia8Decay",pDEBUG) << "   iChan " << iChan << " k=" << k
                                 << " " << dc_daughter[k] << " vs. " << py_daughter[k];
      if (dc_daughter[k] != py_daughter[k]) {
        match = false;  // mismatch
        break; // we're done with this iChan
      }
    }
    // made it through the vector comparison; have a match?
    if (match) {
      LOG("Pythia8Decay",pDEBUG) << "iChan " << iChan << " match " << (match?"yes":"no");
      return iChan;
    }
  } // loop over iChan

#else
  LOG("Pythia8Decay", pFATAL)
    << "calling GENIE/PYTHIA8 decay without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif

  return -1;
}
//___________________________________________________________________________
void Pythia8Decayer2023::PrintDecayChannelInfo(int pdgc) const
{
  LOG("Decay",pDEBUG)
    << "Real implementation for PrintDecayInfo(" << pdgc << ") for "
    << typeid(*this).name();

  TParticlePDG * p = PDGLibrary::Instance()->Find(pdgc);
  if ( ! p ) {
    LOG("Pythia8Decay",pWARN) << "PDGLibrary doesn't know about pdg " << pdgc;
    return;
  }

  // Get the one-and-only instance of Pythia8 that GENIE will use
  Pythia8::Pythia* gPythia = genie::Pythia8Singleton::Instance()->Pythia8();
  if ( ! gPythia->particleData.isParticle(pdgc) ) {
    LOG("Pythia8Decay",pWARN) << "Pythia8 doesn't know about pdg " << pdgc;
    return;
  }

  auto py8_p  = gPythia->particleData.particleDataEntryPtr(pdgc);

  std::ostringstream ss;
  ss << "\n";
  double brsum=0.0, brsum_py8=0.0;
  std::set<int> py8matches;
  for (int ic=0; ic < p->NDecayChannels(); ++ic) {
    TDecayChannel * dch = p->DecayChannel(ic);
    double br     = dch->BranchingRatio();
    double br_py8     =  0;
    int    onMode_py8 = -1;
    // does this PDGLibrary channel map to something in Pythia8?
    int ic_py8 = FindPythiaDecayChannel(pdgc,dch); // -1 if no match
    if (ic_py8>=0) {
      auto py8_dk = py8_p->channel(ic_py8);
      br_py8      = py8_dk.bRatio();
      onMode_py8  = py8_dk.onMode();
    }
    brsum     += br;
    brsum_py8 += br_py8;
    ss << "    " << std::setw(2) << ic
       << " " << std::fixed << std::setprecision(8) << br << " [ ";
    if (ic_py8<0) {
      ss << "-- - ----------";
    } else {
      py8matches.insert(ic_py8);
      ss << std::setw(2) << ic_py8 << " "
         << std::setw(1) << onMode_py8 << " "
         << std::setprecision(8) << br_py8;
    }
    ss << " ] ==> ";
    for (int kdau=0; kdau < dch->NDaughters(); ++kdau) {
      int pdgc_dau = dch->DaughterPdgCode(kdau);
      TParticlePDG * p_dau = PDGLibrary::Instance()->Find(pdgc_dau);
      ss << p_dau->GetName() << "(" << pdgc_dau << ")";
      if ( kdau < dch->NDaughters() -1 ) ss << " + ";
    }
    ss << "\n";
  }
  ss << " BRsum " << std::setprecision(8) << brsum
     << " BRsum_py8 " << std::setprecision(8) << brsum_py8 << "\n";
  ss <<    "The following pythia channels are not accessible from InhibitDecay/Particle:" << "\n";
  //std::set<int> py8missing;
  int npy8chan = py8_p->sizeChannels();
  for (int ic_py8=0; ic_py8 < npy8chan; ++ic_py8) {
    //c++20 if ( ! py8matches.contains(ic_py8) ) {
    if ( py8matches.find(ic_py8) == py8matches.end() ) {
      //py8missing.insert(ic_py8);
      ss << " " << ic_py8;
    }
  }

  LOG("Pythia8Decay",pNOTICE)
    << "\n PDGLibrary decay channels for " <<  p->GetName() << "(" << pdgc << ")"
    << "\n iChan BRatio [ pythia8chan onMode BRatio ] "
    << ss.str() << "\n";

  LOG("Pythia8Decay",pNOTICE)
    << " Pythia8 view of decays for " <<  py8_p->name() << "(" << pdgc << ")";
  gPythia->particleData.list(pdgc);

}
//___________________________________________________________________________
