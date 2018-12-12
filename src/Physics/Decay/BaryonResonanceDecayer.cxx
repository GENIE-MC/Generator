//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TDecayChannel.h>
#include <TMath.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Decay/BaryonResonanceDecayer.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer() :
Decayer("genie::BaryonResonanceDecayer")
{
  this->Initialize();
}
//____________________________________________________________________________
BaryonResonanceDecayer::BaryonResonanceDecayer(string config) :
Decayer("genie::BaryonResonanceDecayer", config)
{
  this->Initialize();
}
//____________________________________________________________________________
BaryonResonanceDecayer::~BaryonResonanceDecayer()
{

}
//____________________________________________________________________________
void BaryonResonanceDecayer::ProcessEventRecord(GHepRecord * event) const
{
  LOG("ResonanceDecay", pINFO)
    << "Running resonance decayer "
    << ((fRunBefHadroTransp) ? "*before*" : "*after*") << " FSI";

  // Loop over particles, find unstable ones and decay them
  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  int ipos = -1;

  while( (p = (GHepParticle *) piter.Next()) ) {

    ipos++;
    LOG("ResonanceDecay", pDEBUG) << "Checking: " << p->Name();

    int pdg_code = p->Pdg();
    GHepStatus_t status_code = p->Status();

    if(!this->IsHandled  (pdg_code)) continue;
    if(!this->ToBeDecayed(pdg_code, status_code)) continue;

    LOG("ResonanceDecay", pNOTICE)
          << "Decaying unstable particle: " << p->Name();

    bool decayed = this->Decay(ipos, event);
    assert(decayed); // handle this more graciously and throw an exception
  } // loop over particles

  LOG("ResonanceDecay", pNOTICE)
     << "Done finding & decaying baryon resonances";
}
//____________________________________________________________________________
bool BaryonResonanceDecayer::Decay(
  int decay_particle_id, GHepRecord * event) const
{
  // Reset previous decay weight
  fWeight = 1.;

  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return false;

  // Select a decay channel
  TDecayChannel * selected_decay_channel =
     this->SelectDecayChannel(decay_particle_id, event);
  if(!selected_decay_channel) return false;

  // Decay the exclusive state and copy daughters in the event record
  this->DecayExclusive(decay_particle_id, event, selected_decay_channel);

  // Update the event weight for each weighted particle decay
  double weight = event->Weight() * fWeight;
  event->SetWeight(weight);

  // Mark input particle as a 'decayed state' & add its daughter links
  decay_particle->SetStatus(kIStDecayedState);

  return true;
}
//____________________________________________________________________________
TDecayChannel * BaryonResonanceDecayer::SelectDecayChannel(
   int decay_particle_id, GHepRecord * event) const
{
  // Get particle to be decayed
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return 0;

  // Get the particle 4-momentum and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Find the particle in the PDG library & quit if it does not exist
  TParticlePDG * mother =
     PDGLibrary::Instance()->Find(decay_particle_pdg_code);
  if(!mother) {
     LOG("ResonanceDecay", pERROR)
        << "\n *** The particle with PDG code = " << decay_particle_pdg_code
         << " was not found in PDGLibrary";
     return 0;
  }
  LOG("ResonanceDecay", pINFO)
    << "Decaying a " << mother->GetName()
    << " with P4 = " << utils::print::P4AsString(&decay_particle_p4);

  // Get the resonance mass W (generally different from the mass associated
  // with the input PDG code, since the it is produced off the mass shell)
  double W = decay_particle_p4.M();
  LOG("ResonanceDecay", pINFO) << "Available mass W = " << W;

  // Get all decay channels
  TObjArray * original_decay_list = mother->DecayList();

  unsigned int nch = original_decay_list -> GetEntries();
  LOG("ResonanceDecay", pINFO)
    << mother->GetName() << " has: " << nch << " decay channels";

  // Loop over the decay channels (dc) and write down the branching
  // ratios to be used for selecting a decay channel.
  // Since a baryon resonance can be created at W < Mres, explicitly
  // check and inhibit decay channels for which W > final-state-mass

  double BR[nch], tot_BR = 0;

  bool is_delta =
    ( decay_particle_pdg_code ==  kPdgP33m1232_Delta0 ||
      decay_particle_pdg_code == -kPdgP33m1232_Delta0 ||
      decay_particle_pdg_code ==  kPdgP33m1232_DeltaP ||
      decay_particle_pdg_code == -kPdgP33m1232_DeltaP );

  TObjArray * actual_decay_list = nullptr ;

  if ( is_delta )
    actual_decay_list = EvolveDeltaBR( decay_particle_pdg_code, original_decay_list, W ) ;
  else
    actual_decay_list = original_decay_list ;

  for(unsigned int ich = 0; ich < nch; ich++) {

    TDecayChannel * ch = (TDecayChannel *) actual_decay_list -> At(ich);

    if ( this->FinalStateMass(ch) < W ) {

      SLOG("ResonanceDecay", pDEBUG)
                << "Using channel: " << ich
                << " with final state mass = " << fsmass << " GeV";

      tot_BR += ch->BranchingRatio();

    } else {
      SLOG("ResonanceDecay", pINFO)
                << "Suppresing channel: " << ich
                << " with final state mass = " << fsmass << " GeV";
    } // final state mass

    BR[ich] = tot_BR;
  }//channel loop

  if ( is_delta )
    delete actual_decay_list ;

  if( tot_BR <= 0. ) {
    SLOG("ResonanceDecay", pWARN)
          << "None of the " << nch << " decay channels is available @ W = " << W;
    return 0;
  }

  // Select a resonance based on the branching ratios
  unsigned int ich = 0, sel_ich; // id of selected decay channel
  RandomGen * rnd = RandomGen::Instance();
  double x = tot_BR * rnd->RndDec().Rndm();
  do {
    sel_ich = ich;
  } while (x > BR[ich++]);

  TDecayChannel * sel_ch = (TDecayChannel *) decay_list->At(sel_ich);

  LOG("ResonanceDecay", pINFO)
    << "Selected " << sel_ch->NDaughters() << "-particle decay channel ("
    << sel_ich << ") has BR = " << sel_ch->BranchingRatio();

  return sel_ch;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::DecayExclusive(
  int decay_particle_id, GHepRecord * event, TDecayChannel * ch) const
{
  // Find the particle to be decayed in the event record
  GHepParticle * decay_particle = event->Particle(decay_particle_id);
  if(!decay_particle) return;

  // Get the decayed particle 4-momentum, 4-position and PDG code
  TLorentzVector decay_particle_p4 = *(decay_particle->P4());
  TLorentzVector decay_particle_x4 = *(decay_particle->X4());
  int decay_particle_pdg_code = decay_particle->Pdg();

  // Get the final state mass spectrum and the particle codes
  // for the selected decay channel
  unsigned int nd = ch->NDaughters();

  int    pdgc[nd];
  double mass[nd];

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {
     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     pdgc[iparticle] = daughter_code;
     mass[iparticle] = daughter->Mass();

     SLOG("ResonanceDecay", pINFO)
         << "+ daughter[" << iparticle << "]: "
         << daughter->GetName() << " (pdg-code = "
         << pdgc[iparticle] << ", mass = " << mass[iparticle] << ")";
  }

  // Check whether the expected channel is Delta->pion+nucleon
  bool is_delta = (decay_particle_pdg_code == kPdgP33m1232_DeltaPP ||
                   decay_particle_pdg_code == kPdgP33m1232_DeltaP  ||
                   decay_particle_pdg_code == kPdgP33m1232_Delta0);

  bool is_delta_N_Pi_decay = is_delta && this->IsPiNDecayChannel(ch);

  // Decay the resonance using an N-body phase space generator
  // The particle will be decayed in its rest frame and then the daughters
  // will be boosted back to the original frame.

  bool is_permitted = fPhaseSpaceGenerator.SetDecay(decay_particle_p4, nd, mass);
  assert(is_permitted);

  // Find the maximum phase space decay weight
  // double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int i=0; i<50; i++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  LOG("ResonanceDecay", pINFO)
    << "Max phase space gen. weight for current decay: " << wmax;

  if(fGenerateWeighted)
  {
    // Generating weighted decays
    // Do a single draw of momentum 4-vectors and then stop,
    // taking into account the weight for this particular draw
    double w = fPhaseSpaceGenerator.Generate();
    fWeight *= TMath::Max(w/wmax, 1.);
  }
  else
  {
    // Generating un-weighted decays
    RandomGen * rnd = RandomGen::Instance();
    wmax *= 2;
    bool accept_decay=false;
    unsigned int itry=0;

    while(!accept_decay)
    {
      itry++;
      assert(itry<kMaxUnweightDecayIterations);

      double w  = fPhaseSpaceGenerator.Generate();
      double gw = wmax * rnd->RndDec().Rndm();

      if(w>wmax) {
         LOG("ResonanceDecay", pWARN)
            << "Current decay weight = " << w << " > wmax = " << wmax;
      }
      LOG("ResonanceDecay", pINFO)
        << "Current decay weight = " << w << " / R = " << gw;

      accept_decay = (gw<=w);

      // Extra logic that applies only for Delta -> N + pi
      if( accept_decay && is_delta_N_Pi_decay ) {

        // We don't want the decay Delta -> Pi + N to be isotropic in the Delta referemce frame
        // as generated by the simple phase space generator.
        // In order to sample pion distribution according to W(Theta) in the Delta resonance decay,
        // we make use of the following.
        // For each event generated from a Delta -> N + Pi event with Pi emitted at
        // at angle Theta (in the Delta rest frame), attach a random number to it.
        // then we calculate W(Theta) according to the definition
        // W(Theta) = 1 − P[ 3/2 ] x L_2(cos Theta) + P[ 1/2 ] x L_2(cos Theta)
        // where
        // L_2 is the second Legendre polynomial L_2(x) = (3x^2 -1)/2
        // and P[3/2] and P[1/2] have to some up to 1.
        // Each possible final state is used to evaluate (Theta),
        // then a random number is thrown, if the the random number is higher than W(Theta) drop this event and go
        // back to re-generate an event and repeat the procedure.
        // Otherwise keep this event to the record.

        // Locate the pion in the decay products
        // at this point we already know that the pion is unique so the first pion we find is our pion
        unsigned int pi_id = 0 ;
        unsigned int nd = ch->NDaughters();

        for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

          if ( genie::pdg::IsPion( ch->DaughterPdgCode(iparticle) ) ) {
            pi_id = iparticle ;
            break ;
          }
        }//iparticle

        TLorentzVector * lab_pion = fPhaseSpaceGenerator.GetDecay(pi_id);
        TLorentzVector pion( lab_pion -> Px(), lab_pion -> Py(), lab_pion -> Pz(), lab_pion -> Energy() ) ;

        pion.Boost(-decay_particle_p4.BoostVector() );  // this gives us the pion in the Delta reference frame

        double costheta = pion.CosTheta() ;
        double legendre_2 = 0.5*(3.*costheta*costheta -1.);

        double w_theta = 1. - fProb32 * legendre_2 + fProb12 * legendre_2 ;

        double aidrnd = 1.25 * rnd->RndDec().Rndm();

        if ( w_theta < aidrnd )
          accept_decay = false ;

      } //if it is a Delta -> N + Pi

    }//accept_decay

  }//fGenerateWeighted

  // A decay was generated - Copy to the event record

  // Check whether the interaction is off a nuclear target or free nucleon
  // Depending on whether this module is run before or after the hadron
  // transport module it would affect the daughters status code
  GHepParticle * target_nucleus = event->TargetNucleus();
  bool in_nucleus = (target_nucleus!=0);

  // Loop over daughter list and add corresponding GHepParticles
  for(unsigned int id = 0; id < nd; id++) {

     int daughter_pdg_code = pdgc[id];

     TLorentzVector * daughter_p4 = fPhaseSpaceGenerator.GetDecay(id);
     LOG("ResonanceDecay", pDEBUG)
        << "Adding daughter particle with PDG code = " << pdgc[id]
        << " and mass = " << mass[id] << " GeV";

     bool is_hadron = pdg::IsHadron(daughter_pdg_code);
     bool hadron_in_nuc = (in_nucleus && is_hadron && fRunBefHadroTransp);

     GHepStatus_t daughter_status_code = (hadron_in_nuc) ?
          kIStHadronInTheNucleus : kIStStableFinalState;

     event->AddParticle(
       daughter_pdg_code, daughter_status_code, decay_particle_id,-1,-1,-1,
       *daughter_p4, decay_particle_x4);
  }

}
//____________________________________________________________________________
double BaryonResonanceDecayer::DealsDeltaNGamma(
  int decay_particle_pdg_code, int ichannel, double W) const
{

/*
 * The branching rations of the Delta in Pions or in N gammas are not constant.
 * They depend on the actual mass of the decaying delta (W) they need to be evolved accordingly.
 * This method tweaks the Delta branching ratios as a function of the W and
 * returns the proper one depending on the specific decay channel.
 */


  if (W <= genie::constants::kNucleonMass+genie::genie::constants::kPi0Mass) {

    if (ichannel == 0) {return 0;} // ichannel = 0,1,2 has to match
                                   // the channel order in genie_pdg_table.dat
    if (ichannel == 1) {return 0;}
    if (ichannel == 2) {return 1;}
  }
  else
  {
    // LIBO: Add detailed explanation for the calculation you are doing here
    //
    //
    //
    //
    //
    //

    // Getting the delta resonance from GENIE database
    Resonance_t res = genie::utils::res::FromPdgCode( decay_particle_pdg_code ) ;

    // get the width of the delta and obtain the width of the decay in Pi+N and gamma+N
    // evaluated at the nominal mass of the delta
    double width0 = genie::utils::res::Width( res ) ;
    double widPi0   = width0*fPionBR;
    double widgamma0= width0*fGammaBR;

    // Possible optimization here: we could use the actual masses of the final state particles
    double m      = genie::utils::res::Mass( res ) ;
    double m_2   = TMath::Power(m,      2);
    double mN_2  = TMath::Power(genie::constants::kNucleonMass,     2);
    double W_2   = TMath::Power(W,      2);
    double m_aux1= TMath::Power(genie::constants::kNucleonMass+genie::genie::constants::kPi0Mass, 2);
    double m_aux2= TMath::Power(genie::constants::kNucleonMass-genie::genie::constants::kPi0Mass, 2);

    double pPiW    = TMath::Sqrt((W_2-m_aux1)*(W_2-m_aux2))/(2*W);
    double pPim    = TMath::Sqrt((m_2-m_aux1)*(m_2-m_aux2))/(2*m);
    double EgammaW = (W_2-mN_2)/(2*W);
    double Egammam = (m_2-mN_2)/(2*m);
    double TPiW    = TMath::Power(pPiW, 3);
    double TPim    = TMath::Power(pPim, 3);
    double fgammaW = 1./(TMath::Power(1+EgammaW*EgammaW/fFFScaling, 2));
    double fgammam = 1./(TMath::Power(1+Egammam*Egammam/fFFScaling, 2));

    double Rinverse =
        widPi0*TMath::Power(Egammam, 3)*TMath::Power(fgammam, 2)*TPiW /
	      (widgamma0*TMath::Power(EgammaW, 3)*TMath::Power(fgammaW, 2)*TPim);
    double BRPi = Rinverse/(1+Rinverse);
    double BRgamma = 1/(1+Rinverse);

    double BRPi01   = 0.667002; // LIBO: ditto
    double BRPi02   = 0.332998; // LIBO: ditto

    // Delta0 or Delta0_bar
    if (decay_particle_pdg_code ==  kPdgP33m1232_Delta0 ||
        decay_particle_pdg_code == -kPdgP33m1232_Delta0)
    {
   	   if (ichannel == 0) { return BRPi*BRPi02; }
	     if (ichannel == 1) { return BRPi*BRPi01; }
	     if (ichannel == 2) { return BRgamma;     }
    }
    // Delta+ or anti_Delta+
    else
    if (decay_particle_pdg_code ==  kPdgP33m1232_DeltaP ||
        decay_particle_pdg_code == -kPdgP33m1232_DeltaP)
    {
  	   if (ichannel == 0) { return BRPi*BRPi01; }
	     if (ichannel == 1) { return BRPi*BRPi02; }
	     if (ichannel == 2) { return BRgamma;     }
    }
    else
    {
      LOG("ResonanceDecay", pWARN)
         << "Mother particle (PDG code = " << decay_particle_pdg_code
         << ") is not Delta+ or Delta0!";
  	}
  }//W
  return 0;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::EvolveDeltaBR(int dec_part_pdgc, TDecayChannel * ch, double W) const {

  /*
   * The branching rations of the Delta in Pions or in N gammas are not constant.
   * They depend on the actual mass of the decaying delta (W) they need to be evolved accordingly.
   * This method tweaks the Delta branching ratios as a function of the W and
   * returns the proper one depending on the specific decay channel.
   */

  // identify the decay channel
  // The delta decays only in 3 ways
  // Delta -> Charged Pi + N
  // Delta -> Pi0 + N
  // Delta -> Gamma + N

  // They have evolution as a function of W that are different if the final state has pions or not
  // so having tagged the pion is enough for the purpose of this method.

  bool has_pion = false ;
  int pion_id = -1 ;
  int nucleon_id = -1 ;
  unsigned int nd = ch -> NDaughters() ;
  for( int i = 0 ; i < nd; ++i ) {
    if ( genie::pdg::IsPion( ch -> DaughterPdgCode(i) ) ) {
      has_pion = true ;
      pion_id = i ;
    }

    if ( genie::pdg::IsNucleon( ch -> DaughterPdgCode(i) ) ) {
      nucleon_id = i ;
    }
  }


  // The first and most trivial evolution of the BR as a function of W
  // is that if W is lower then mass of Pi0 and the corresponding nucleon
  // The only possible decay is the one with gamma, therefore the BR goes
  // to either 0 or 1 depending if there is pion or not in the final state.
  // The reason why we check only with pi0 is because reguardless of the
  // charge of the Delta, the ligther final state is always the one with the pi0
  // For simplicity the check is done on the generic nucleon mass

  if (W < genie::constants::kPi0Mass + genie::constants::kNucleonMass ) {

    if ( has_pion ) return 0. ;
    else return 1. ;

   }

  // there is a region of a few MeV in which the Pi0 decay is allowed by the W and the other is not
  // this condition is treated as both the decays were available

  // At this point, W is high enough to assume the decay of the delta in both N+pi or N+gamma
  // This requires the amplitude of both decays to be scaled according to W
  // The amplitude dependencies of W scales with the momentum of the pion or the photon respectivelly
  // following these relationships
  //
  //                              (p_pi(W))^3
  //  Ampl_pi(W) = Ampl_pi(def)x---------------
  //                             (p_pi(def))^3
  //
  //
  //                              (p_ga(W))^3       (F_ga(W))^2
  //  Ampl_ga(W) = Ampl_ga(def)x--------------- x ---------------
  //                             (p_ga(def))^3     (F_ga(def))^2
  //
  // where the "def" stand for the nominal value of the Delta mass.
  //  - pi_* are the momentum of the gamma and of the pion coming from the decay
  //  - F_ga is the form factor
  //
  // So the new amplitudes are evaluated and the proper value is returned

  // Getting the delta resonance from GENIE database
   Resonance_t res = genie::utils::res::FromPdgCode( dec_part_pdgc ) ;

   // get the width of the delta and obtain the width of the decay in Pi+N and gamma+N
   // evaluated at the nominal mass of the delta
   double defWidth   = genie::utils::res::Width( res ) ;

   double m = genie::utils::res::Mass( res ) ;
   double m_2   = TMath::Power(m, 2);

   double mN = genie::pdg::IsProton( ch -> DaughterPdgCode( nucleon_id ) ) ?  genie::constants::kProtonMass : genie::constants::kNucleonMass ;
   double mN_2  = TMath::Power( mN,     2);

   double mPion = TMath::Abs( ch -> DaughterPdgCode( pion_id ) ) == kPdgPiP ? genie::constants::kPionMass : genie::constants::kPi0Mass ;
   double m_aux1= TMath::Power( mN + mPion, 2) ;
   double m_aux2= TMath::Power( mN - mPion, 2) ;

   double W_2   = TMath::Power(W,      2);

   // momentum of the pion in the Delta reference frame
   double pPi_W    = TMath::Sqrt((W_2-m_aux1)*(W_2-m_aux2))/(2*W);  // at W
   double pPi_m    = TMath::Sqrt((m_2-m_aux1)*(m_2-m_aux2))/(2*m);  // at the default Delta mass

   double TPi_W    = TMath::Power(pPi_W, 3);
   double TPi_m    = TMath::Power(pPi_m, 3);

   // momentum of the photon in the Delta Reference frame = Energy of the photon
   double Egamma_W = (W_2-mN_2)/(2*W);  // at W
   double Egamma_m = (m_2-mN_2)/(2*m);  // at the default Delta mass

   // form factor of the photon production
   double fgamma_W = 1./(TMath::Power(1+Egamma_W*Egamma_W/fFFScaling, 2));
   double fgamma_m = 1./(TMath::Power(1+Egamma_m*Egamma_m/fFFScaling, 2));


   double defPiWidth = width0*fPionBR;
   double defGaWidth = width0*fGammaBR;


}
//____________________________________________________________________________

double BaryonResonanceDecayer::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
double BaryonResonanceDecayer::FinalStateMass(TDecayChannel * ch) const
{
// Computes the total mass of the final state system

  double mass = 0;
  unsigned int nd = ch->NDaughters();

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

     int daughter_code = ch->DaughterPdgCode(iparticle);
     TParticlePDG * daughter = PDGLibrary::Instance()->Find(daughter_code);
     assert(daughter);

     double md = daughter->Mass();

     // hack to switch off channels giving rare  occurences of |1114| that has
     // no decay channels in the pdg table (08/2007)
     if(TMath::Abs(daughter_code) == 1114) {
         LOG("ResonanceDecay", pNOTICE)
            << "Disabling decay channel containing resonance 1114";;
         md = 999999999;
     }
     mass += md;
  }
  return mass;
}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsPiNDecayChannel(TDecayChannel * ch) const
{
  if(!ch) return false;

  unsigned int nd = ch->NDaughters();
  if(nd != 2) return false;

  int  npi    = 0;
  int  nnucl  = 0;

  for(unsigned int iparticle = 0; iparticle < nd; iparticle++) {

    int daughter_code = ch->DaughterPdgCode(iparticle);

    if( genie::pdg::IsPion( daughter_code ) )
      npi++;

    if ( genie::pdg::IsNucleon(daughter_code ) )
      nnucl++;

  }//iparticle

  if(npi == 1 && nnucl == 1) return true;

  return false;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::Initialize(void) const
{

}
//____________________________________________________________________________
bool BaryonResonanceDecayer::IsHandled(int pdg_code) const
{
  bool is_handled = utils::res::IsBaryonResonance(pdg_code);

  LOG("ResonanceDecay", pDEBUG)
      << "Can decay particle with PDG code = " << pdg_code
      << "? " << ((is_handled)? "Yes" : "No");

  return pdg_code;
}
//____________________________________________________________________________
void BaryonResonanceDecayer::InhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
void BaryonResonanceDecayer::UnInhibitDecay(int pdgc, TDecayChannel * dc) const
{
  if(! this->IsHandled(pdgc)) return;
  if(!dc) return;

  //
  // Not implemented
  //
}
//____________________________________________________________________________
void BaryonResonanceDecayer::LoadConfig(void) {

  Decayer::LoadConfig() ;

  this->GetParam( "Prob32", fProb32 ) ;
  fProb12 = 1. - fProb32 ;

  this->GetParam( "TotGammaBR", fGammaBR ) ;
  fPionBR = 1. - fGammaBR ;

  fWidthPi_0 =    delta_width * TotNPiBR ;
  fWidthGamma_0 = delta_width * TotGammaBR ;

  this -> GetParam( "OnePiBR", fBRPi01 ) ;
  fBRPi02 = 1. - fBRPi01 ;

  fDeltaMass2   = TMath::Power( fDeltaMass, 2) ;
  fNucleonMass2 = TMath::Power( genie::constants::kNucleonMass, 2) ;
  fMAux1 = TMath::Power( genie::constants::kNucleonMass + genie::constants::kPi0Mass , 2) ;
  fMAux2 = TMath::Power( genie::constants::kNucleonMass - genie::constants::kPi0Mass , 2) ;

  this -> GetParam( "FFScaling", fFFScaling ) ;


}

