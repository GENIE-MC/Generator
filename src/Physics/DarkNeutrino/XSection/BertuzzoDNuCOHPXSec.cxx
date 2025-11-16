//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
  University of Sussex

  Costas Andreopoulos <c.andreopoulos \at cern.ch>
  University of Liverpool
*/
//____________________________________________________________________________


#include <TMath.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/DarkNeutrino/XSection/BertuzzoDNuCOHPXSec.h"
#include "Physics/DarkNeutrino/XSection/EngelFormFactor.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
BertuzzoDNuCOHPXSec::BertuzzoDNuCOHPXSec() :
XSecAlgorithmI("genie::BertuzzoDNuCOHPXSec")
{

}
//____________________________________________________________________________
BertuzzoDNuCOHPXSec::BertuzzoDNuCOHPXSec(string config) :
XSecAlgorithmI("genie::BertuzzoDNuCOHPXSec", config)
{

}
//____________________________________________________________________________
BertuzzoDNuCOHPXSec::~BertuzzoDNuCOHPXSec()
{

}
//____________________________________________________________________________
double BertuzzoDNuCOHPXSec::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  // User inputs to the calculation
  const int nu_pdg = init_state.ProbePdg();
  const double E  = init_state.ProbeE(kRfLab); // neutrino energy, units: GeV
  const double Q2 = kinematics.Q2(); // momentum transfer, units: GeV^2
  const double TE =  kinematics.HadSystP4().E(); // energy of the target
  const unsigned int Z = target.Z(); // number of protons

  // select the mixing depending on the incoming neutrino
  unsigned short nu_i = 3;
  if( pdg::IsNuE( TMath::Abs( nu_pdg ) ) ) nu_i = 0;
  else if ( pdg::IsNuMu( TMath::Abs( nu_pdg ) ) ) nu_i = 1;
  else if ( pdg::IsNuTau( TMath::Abs( nu_pdg ) ) ) nu_i = 2;
  const double DTheta2 = fMixing2s[nu_i] * fMixing2s[3] ;

  // Target atomic mass number and mass calculated from inputs
  const double M = target.Mass(); // units: GeV

  const double FF = fFF->FormFactor(Q2, target);
  const double TT = TE - M;

  // auxiliary variables
  const double E2  = E * E;
  const double Z2 = Z * Z;
  const double FF2 = FF * FF;
  const double TTDiff = TT - M;

  const double const_factor = 2* constants::kPi * constants::kAem ;
  const double model_params = fEps2 * DTheta2 * fAlpha_D ;

  const double num_fact1 = FF2 * Z2;
  const double num_fact21 = fDNuMass2 * (TTDiff - 2.*E);
  const double num_fact22 = 2. * M * (2.*E2 - 2.*TT*E + TT*TTDiff);
  const double den_fact1 = 1. / (E2);
  const double den_fact2 = TMath::Power((fDMediatorMass2 + 2.*TT*M), -2.);

  if(kps == kPSEDNufE) {
    const double xsec = const_factor * model_params * num_fact1 *
      (num_fact21 + num_fact22) * den_fact1 * den_fact2;

    return xsec;
  }
  return 0.;
}
//____________________________________________________________________________
double BertuzzoDNuCOHPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BertuzzoDNuCOHPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();
  if ( ! proc_info.IsCoherentElastic() ) return false;
  if ( ! proc_info.IsDarkNeutralCurrent() ) return false ;

  const InitialState & init_state = interaction->InitState();
  if ( ! pdg::IsNeutrino( TMath::Abs( init_state.ProbePdg() ) ) ) return false ;

  const Target & target = init_state.Tgt();
  if( ! target.IsNucleus() ) return false ;

  return true;
}
//____________________________________________________________________________
bool BertuzzoDNuCOHPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  if(!interaction->PhaseSpace().IsAboveThreshold()) return false;

  const double E  = interaction->InitState().ProbeE(kRfLab);
  const double M = interaction->InitState().Tgt().Mass();
  const TLorentzVector& DNu = interaction->Kine().FSLeptonP4();

  const double tl = DNu.E()*(M+E) - E*M - 0.5*fDNuMass2;
  const double tr = E * DNu.P();

  if(tl < -1.*tr) return false;
  if(tl >  tr) return false;

  return true;

}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BertuzzoDNuCOHPXSec::LoadConfig(void)
{

  bool good_configuration = true ;

  double DKineticMixing = 0.;
  this->GetParam("Dark-KineticMixing", DKineticMixing);
  fEps2 = DKineticMixing * DKineticMixing;

  bool force_unitarity = false ;
  GetParam( "Dark-Mixing-ForceUnitarity", force_unitarity ) ;

  unsigned int n_min_mixing = force_unitarity ? 3 : 4 ;

  std::vector<double> DMixing2s;  // |U_{\alpha 4}|^2
  this->GetParamVect("Dark-Mixings2", DMixing2s);

  // check whether we have enough mixing elements
  if ( DMixing2s.size () < n_min_mixing ) {
    good_configuration = false ;
    LOG("BertuzzoDNuCOH", pERROR )
      << "Not enough mixing elements specified, only specified "
      << DMixing2s.size() << " / " << n_min_mixing ;
  }

  double tot_mix = 0.;
  for( unsigned int i = 0; i < n_min_mixing ; ++i ) {
    if ( DMixing2s[i] < 0. ) {
      good_configuration = false ;
      LOG("BertuzzoDNuCOH", pERROR )
        << "Mixing " << i << " non positive: " << DMixing2s[i] ;
      continue ;
    }
    tot_mix += fMixing2s[i] = DMixing2s[i] ;
  }

  if ( force_unitarity ) {
    fMixing2s[3] = 1. - tot_mix ;
  }
  if ( DMixing2s[3] < 0. ) {
    good_configuration = false ;
    LOG("BertuzzoDNuCOH", pERROR )
      << "Mixing D4 non positive: " << DMixing2s[3] ;
  }
  
  this->GetParam("Dark-Alpha", fAlpha_D);

  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);
  fDNuMass2 = fDNuMass * fDNuMass;

  fDMediatorMass = 0.;
  this->GetParam("Dark-MediatorMass", fDMediatorMass);
  fDMediatorMass2 = fDMediatorMass * fDMediatorMass;

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
  fFF = dynamic_cast<const EngelFormFactor *> (this->SubAlg("FormFactor"));
  assert(fFF);

  if ( ! good_configuration ) {
    LOG("BertuzzoDNuCOH", pFATAL ) << "Wrong configuration. Exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
