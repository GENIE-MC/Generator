//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
 University of Sussex

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cmath>
#include <numeric>

#include <TMath.h>

#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Decay/DarkSectorDecayer.h"
#include "Math/GSLMinimizer.h"
#include <Math/WrappedParamFunction.h>

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;

//____________________________________________________________________________
DarkSectorDecayer::DarkSectorDecayer() :
  EventRecordVisitorI("genie::DarkSectorDecayer")
{

}
//____________________________________________________________________________
DarkSectorDecayer::DarkSectorDecayer(string config) :
  EventRecordVisitorI("genie::DarkSectorDecayer", config)
{

}
//____________________________________________________________________________
DarkSectorDecayer::~DarkSectorDecayer()
{

}
//____________________________________________________________________________
void DarkSectorDecayer::ProcessEventRecord(GHepRecord * event) const
{
  LOG("DarkSectorDecayer", pINFO)
    << "Running dark sector decayer ";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {
    ipos++;
    LOG("DarkSectorDecayer", pDEBUG) << "Checking: " << p->Name();

    if(!this->ToBeDecayed(*p)) continue;

    GHepParticle&  mother = *p; // change the name now we know it will decay
    std::vector<DarkSectorDecayer::DecayChannel> dcs;
    int pdg_code = mother.Pdg();
    if(pdg_code == kPdgDNuMediator){
      dcs = DarkMediatorDecayChannels();
    }
    else if(pdg_code == kPdgDarkNeutrino ||
            pdg_code == kPdgAntiDarkNeutrino){
      dcs = DarkNeutrinoDecayChannels(pdg_code);
    }

    // for ( const auto & dc : dcs ) {

    // std::cout << "Decay amplitude: " << dc.second << " GeV "
    //           << " -> " << 1. / dc.second / units::second << "s for channel [" ;
    //   for ( const auto & p : dc.first ) {
    //       std::cout << p << "  " ;
    //   }
    //   std::cout << "]" << std::endl ;
    // }

    double total_amplitude = std::accumulate(dcs.begin(), dcs.end(), 0.,
                                             [](double total,
                                                const DarkSectorDecayer::DecayChannel& dc)
                                               {return total + dc.second;});

    int dcid = SelectDecayChannel(dcs, total_amplitude);
    std::vector<GHepParticle> daughters = Decay(mother, dcs[dcid].first);
    SetSpaceTime(daughters, mother, total_amplitude);

    for(auto & daughter: daughters){
      daughter.SetFirstMother(ipos);
      event->AddParticle(daughter);
    }
  }

  LOG("DarkSectorDecayer", pNOTICE)
     << "Done finding & decaying dark sector particles";

}
//____________________________________________________________________________
std::vector<GHepParticle> DarkSectorDecayer::Decay(
  const GHepParticle & mother,
  const std::vector<int> & pdg_daughters) const
{
  TLorentzVector mother_p4 = *(mother.P4());
  LOG("DarkSectorDecayer", pINFO)
    << "Decaying a " << mother.GetName()
    << " with P4 = " << utils::print::P4AsString(&mother_p4);

  unsigned int nd = pdg_daughters.size();
  double mass[nd] = {0.};

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {
    TParticlePDG * daughter = PDGLibrary::Instance()->Find(pdg_daughters[iparticle]);
    assert(daughter);

    mass[iparticle] = daughter->Mass();

    SLOG("DarkSectorDecayer", pINFO)
      << "+ daughter[" << iparticle << "]: "
      << daughter->GetName() << " (pdg-code = "
      << pdg_daughters[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(mother_p4, nd, mass);
  assert(is_permitted);

  // Find the maximum phase space decay weight
  double wmax = -1;
  for(int i=0; i<50; i++) {
    double w = fPhaseSpaceGenerator.Generate();
    wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("DarkSectorDecayer", pINFO)
    << "Max phase space gen. weight for current decay: " << wmax;

  // Generating un-weighted decays
  RandomGen * rnd = RandomGen::Instance();
  wmax *= 2;
  bool accept_decay=false;
  unsigned int itry=0;

  while(!accept_decay){
    itry++;
    assert(itry<kMaxUnweightDecayIterations);

    double w  = fPhaseSpaceGenerator.Generate();
    double gw = wmax * rnd->RndDec().Rndm();

    if(w>wmax) {
      LOG("DarkSectorDecayer", pWARN)
        << "Current decay weight = " << w << " > wmax = " << wmax;
    }
    LOG("DarkSectorDecayer", pINFO)
      << "Current decay weight = " << w << " / R = " << gw;

    accept_decay = (gw<=w);
  }

  // A decay was generated - Copy to the event record
  std::vector<GHepParticle> particles;
  // Loop over daughter list and add corresponding GHepParticles
  for(unsigned int id = 0; id < nd; id++) {
    TLorentzVector * daughter_p4 = fPhaseSpaceGenerator.GetDecay(id);
    LOG("DarkSectorDecayer", pDEBUG)
      << "Adding daughter particle with PDG code = " << pdg_daughters[id]
      << " and mass = " << mass[id] << " GeV";
    GHepStatus_t daughter_status_code = (pdg_daughters[id]==kPdgDNuMediator)
      ? kIStDecayedState : kIStStableFinalState;
    particles.push_back(GHepParticle(pdg_daughters[id], daughter_status_code,
                                     -1, -1, -1, -1,
                                     *daughter_p4, TLorentzVector()));
  }

  return particles;
}
//____________________________________________________________________________
int DarkSectorDecayer::SelectDecayChannel(
  const std::vector<DarkSectorDecayer::DecayChannel> & dcs,
  const double total_amplitude) const
{
  // Select a decay based on the amplitudes
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = total_amplitude * rnd->RndDec().Rndm();
  double partial_sum = 0. ;
  do {
    sel_ich = ich;
    partial_sum += dcs.at(ich++).second;
  } while (x > partial_sum );
  return sel_ich;
}
//____________________________________________________________________________
std::vector<DarkSectorDecayer::DecayChannel> DarkSectorDecayer::DarkMediatorDecayChannels(void) const
{
  // eq (4) and (5) and maybe some other higher order variations
  // TODO DNu: what alpha_D is?

  const double alpha_D = 0.25; // value on the paper
  std::array<int, 3> neutrinos = {kPdgNuE, kPdgNuMu, kPdgNuTau};
  std::array<int, 3> antineutrinos = {kPdgAntiNuE, kPdgAntiNuMu, kPdgAntiNuTau};
  std::vector<DarkSectorDecayer::DecayChannel> dcs;

  for(size_t i=0; i<neutrinos.size(); ++i){
    for(size_t j=0; j<antineutrinos.size(); ++j){// for antineutrinos
      const double decay_width = alpha_D/3. * fMixing2s[i]*fMixing2s[j] * fDMediatorMass;
      dcs.push_back(DecayChannel{{neutrinos[i], antineutrinos[j]}, decay_width});
    }
  }

  if(fDMediatorMass > 2.*PDGLibrary::Instance()->Find(kPdgElectron)->Mass()){
    const double decay_width = kAem*fEps2/3. * fDMediatorMass;
    dcs.push_back(DecayChannel{{kPdgElectron, kPdgPositron}, decay_width});
  }
  // In the future for the decay to muons
  // if(fDMediatorMass > 2.*PDGLibrary::Instance()->Find(kPdgMuon)->Mass()){
  //   const double decay_width = kAem*epsilon2/3. * fDMediatorMass;
  //   dcs.push_back(DecayChannel{{kPdgMuon, kPdgAntiMuon}, decay_width});
  // }
  return dcs;
}
//____________________________________________________________________________
std::vector<DarkSectorDecayer::DecayChannel> DarkSectorDecayer::DarkNeutrinoDecayChannels(int mother_pdg) const
{
  // eq (3) and higher order variations

  const double alpha_D = 0.25; // value on the paper
  std::array<int, 3> neutrinos = {kPdgNuE, kPdgNuMu, kPdgNuTau};
  std::array<int, 3> antineutrinos = {kPdgAntiNuE, kPdgAntiNuMu, kPdgAntiNuTau};
  std::vector<DarkSectorDecayer::DecayChannel> dcs;

  if(fDNuMass > fDMediatorMass){
    for(size_t i=0; i<neutrinos.size(); ++i){
      const double mass2ratio = fDMediatorMass2/fDNuMass2;
      const double p0 = 0.5*alpha_D * fMixing2s[3] * fMixing2s[i];
      const double p1 = fDNuMass*fDNuMass2/fDMediatorMass2;
      const double p2 = 1 - mass2ratio;
      const double p3 = 1 + mass2ratio - 2*mass2ratio*mass2ratio;
      const double decay_width = p0 * p1 * p2 * p3;

      if(mother_pdg == kPdgDarkNeutrino){
        dcs.push_back(DecayChannel{{neutrinos[i], kPdgDNuMediator}, decay_width});
      }
      else if(mother_pdg == kPdgAntiDarkNeutrino){
        dcs.push_back(DecayChannel{{antineutrinos[i], kPdgDNuMediator}, decay_width});
      }
    }
  }
  return dcs;
}
//____________________________________________________________________________
void DarkSectorDecayer::SetSpaceTime(
  std::vector<GHepParticle> & pp,
  const GHepParticle & mother,
  double total_amplitude) const
{
  // TODO DNU: pay attention to units being used in GENIE!
  // convert decay amplitude into time
  const double lifetime =  1e24/units::second/total_amplitude; // units of 10ˆ-24 s

  RandomGen * rnd = RandomGen::Instance();
  double t = rnd->RndDec().Exp(lifetime);

  // t is the decay time in the mother reference frame
  // it needs to be boosted by a factor gamma
  t *= mother.P4() -> Gamma() ;

  // get beta of decaying particle
  const TLorentzVector & mother_X4 = *(mother.X4());
  TVector3 mother_boost = mother.P4()->BoostVector();

  // transport decay_particle with respect to their mother
  double speed_of_light = units::second/units::meter; // this gives us the speed of light in m/s
  TVector3 daughter_position = mother_X4.Vect() + mother_boost * (speed_of_light * t * 1e-9);// in fm

  for(auto & p : pp){
    p.SetPosition(TLorentzVector(daughter_position, mother_X4.T() + t));
  }
}
//____________________________________________________________________________
bool DarkSectorDecayer::ToBeDecayed(const GHepParticle & p) const
{
  GHepStatus_t status_code = p.Status();
  if(status_code != kIStDecayedState) return false;

  int pdg_code = p.Pdg();
  bool is_handled = false;
  if(pdg_code == kPdgDNuMediator ||
     pdg_code == kPdgDarkNeutrino ||
     pdg_code == kPdgAntiDarkNeutrino){
    is_handled = true;
  }

  LOG("DarkSectorDecayer", pDEBUG)
      << "Can decay particle with PDG code = " << pdg_code
      << "? " << ((is_handled)? "Yes" : "No");

  // // Find the particle in the PDG library & quit if it does not exist
  // TParticlePDG * mother = PDGLibrary::Instance()->Find(pdg_code);
  // if(!mother && is_handled) {
  //   LOG("DarkSectorDecayer", pERROR)
  //     << "\n *** The particle with PDG code = " << pdg_code
  //     << " was not found in PDGLibrary";
  //   //exit;
  // }

  return is_handled;
}
//____________________________________________________________________________
void DarkSectorDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DarkSectorDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DarkSectorDecayer::LoadConfig(void)
{

  bool good_configuration = true ;

  double DKineticMixing = 0.;    // \varepsilon
  this->GetParam("Dark-KineticMixing", DKineticMixing);
  fEps2 = DKineticMixing * DKineticMixing;

  bool force_unitarity = false ;
  GetParam( "Dark-Mixing-ForceUnitarity", force_unitarity ) ;

  unsigned int n_min_mixing = force_unitarity ? 3 : 4 ;

  std::vector<double> DMixing2s;  // |U_{\alpha 4}|^2
  this->GetParamVect("Dark-Mixings2", DMixing2s);

  // check whther we go enough mixing elements
  if ( DMixing2s.size () < n_min_mixing ) {

    good_configuration = false ;
    LOG("DarkSectorDecayer", pERROR ) << "Not enough mixing elements specified, only specified "
                                      << DMixing2s.size() << " / " << n_min_mixing ;
  }

  double tot_mix = 0. ;
  for( unsigned int i = 0; i < n_min_mixing ; ++i ) {
    if ( DMixing2s[i] < 0. ) {
      good_configuration = false ;
      LOG("DarkSectorDecayer", pERROR ) << "Mixign " << i << " non positive: " << DMixing2s[i] ;
      continue ;
    }
    tot_mix += fMixing2s[i] = DMixing2s[i] ;
  }

  if ( force_unitarity ) {
    fMixing2s[3] = 1. - tot_mix ;
  }

  double DGaugeCoupling = 0.;   // g_D
  this->GetParam("Dark-GaugeCoupling", DGaugeCoupling);
  fgD2 = DGaugeCoupling * DGaugeCoupling;

  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);
  fDNuMass2 = fDNuMass * fDNuMass;

  fDMediatorMass = 0.;
  this->GetParam("Dark-MediatorMass", fDMediatorMass);
  fDMediatorMass2 = fDMediatorMass * fDMediatorMass;

  // the model is build on the assumption that the mass of the mediator
  // is smaller than the mass of the dark neutrino.
  // For the cross section, this is not a problem
  // The decayer though is sensitive to this as the only known decay amplitude
  // of the dark neutrino requires the mediator in the final state.
  // Until the decay aimplitude in neutrino is not available
  // we need to check that the mass hierarchy is respected
  if ( fDMediatorMass >= fDNuMass ) {
    good_configuration = false ;
    LOG("DarkSectorDecayer", pERROR ) << "Dark mediator mass (" <<  fDMediatorMass
                                      << " GeV) too heavy for the dark neutrino ("
                                      << fDNuMass << " GeV) to decay" ;
  }

  if ( ! good_configuration ) {
    LOG("DarkSectorDecayer", pFATAL ) << "Wrong configuration. Exiting" ;
    exit ( 78 ) ;
  }

}
