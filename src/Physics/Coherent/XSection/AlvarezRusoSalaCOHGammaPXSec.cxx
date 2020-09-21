//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  Jon Sensenig 
  University of Pennsylvania
  July 2020
  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Physics/Coherent/XSection/AlvarezRusoSalaCOHGammaPXSec.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"



using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
AlvarezRusoSalaCOHGammaPXSec::AlvarezRusoSalaCOHGammaPXSec() :
  XSecAlgorithmI("genie::AlvarezRusoSalaCOHGammaPXSec"),
  fFormFactors( nullptr ),
  fXSecIntegrator( nullptr )
{ ; }

//____________________________________________________________________________
AlvarezRusoSalaCOHGammaPXSec::AlvarezRusoSalaCOHGammaPXSec(string config) :
  XSecAlgorithmI("genie::AlvarezRusoSalaCOHGammaPXSec", config),
  fFormFactors( nullptr ),
  fXSecIntegrator( nullptr )
{ ; }
//____________________________________________________________________________
AlvarezRusoSalaCOHGammaPXSec::~AlvarezRusoSalaCOHGammaPXSec()
{ ; }
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::XSec( const Interaction * interaction,
                                           KinePhaseSpace_t kps ) const {

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  if ( kps != kPSEgTlOgfE     && 
       kps != kPSEgTlTgPgfE   && 
       kps != kPSEgOlOgfE ) {
    
    LOG("AlvarezRusoSalaCOHGammaPXSec", pERROR ) 
      << "Requested cross section with invalid phase space: " << KinePhaseSpace::AsString( kps ) ;
    return 0. ;
  } ;

  // the form factor in this interaction is limited to a very narrow region
  // the most efficient way to evaluate the cross section is to first evaluate 
  // the form factor(s) and only then, if they are not zero, proceed to the evaluation
  // of the hadronic currents

  const Kinematics &   kinematics = interaction -> Kine();
  double transferred_Q = sqrt( kinematics.t() ) ;
  int pdg = interaction -> InitState().Tgt().Pdg() ;
  double proton_ff = fFormFactors -> ProtonFF( transferred_Q, pdg ) ;
  double neutron_ff = fFormFactors -> NeutronFF( transferred_Q, pdg ) ;
  
  if ( proton_ff == 0. && neutron_ff == 0. ) return 0. ; 

  utils::math::GTraceContraction H( TotalHadronicCurrent( *interaction, proton_ff, neutron_ff ) ) ;

  const InitialState & init_state = interaction -> InitState();
  
  int nu = init_state.ProbePdg();
  double contraction = pdg::IsNeutrino(nu) ? 
    NeutrinoHadronContraction( *interaction, H ) :
    AntiNeutrinoHadronContraction( *interaction, H ) ;    
  
  double pre_factor = ( 4.*constants::kPi*constants::kAem * constants::kGF2 ) / ( 16*pow(2.*constants::kPi, 5) );

  TLorentzVector gamma_p4 = kinematics.HadSystP4() ; // FS gamma momentum
  TLorentzVector out_lep = kinematics.FSLeptonP4() ; 
   
  double diff_cross_section = pre_factor*contraction*( gamma_p4.E() * out_lep.E() / init_state.ProbeE(kRfLab) ); 

  if ( kps == kPSEgOlOgfE ) return diff_cross_section ;

  // the system is clearly invariant under z rotations -> phi_lepton translation.
  // so that integral can be performed simply multiplying by 2 pi
  diff_cross_section *= 2.*constants::kPi ; 

  // the remaining possible phase spaces are all functions of Theta_l instead of Cos Theta_l
  // so the first jacobian correction can be applied here 
  // teh Jacobian simply being sin_Theta_l
  
  
  diff_cross_section *= sin( out_lep.Theta() ) ;
  
  if ( kps == kPSEgTlOgfE ) 
    return diff_cross_section ;


  // if kps == kPSEgTlTgPgfE 
  diff_cross_section *= sin( gamma_p4.Theta() ) ;
  return diff_cross_section ;   

}
//____________________________________________________________________________
utils::math::GTrace 
AlvarezRusoSalaCOHGammaPXSec::TotalHadronicCurrent( const Interaction & interaction, 
						    double proton_ff, double neutron_ff ) const {

  // sum R's from all currents here
  utils::math::GTrace R_total;
  for ( unsigned int i = 0; i < fCurrents.size(); i++ ) {
    R_total += fCurrents[i] -> R(interaction, proton_ff, neutron_ff );
  }

  return R_total;
}
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::Integral( const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlvarezRusoSalaCOHGammaPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  int nu = init_state.ProbePdg();
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  const ProcessInfo & proc_info = interaction->ProcInfo();

  if ( ! proc_info.IsCoherentProduction() )  return false;
  if ( ! proc_info.IsWeakNC() )              return false;

  const Target & target = init_state.Tgt();
  if ( target.A() <= 1 ) return false;
  if ( ! fFormFactors -> HasNucleus( target.Pdg() ) ) return false ;

  const XclsTag & xcls = interaction -> ExclTag() ;
  if ( xcls.NSingleGammas() != 1 ) return false ;

  return true;
}
//____________________________________________________________________________
bool AlvarezRusoSalaCOHGammaPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const KPhaseSpace& kps = interaction->PhaseSpace();

  if ( ! kps.IsAllowed() ) {
    LOG("AlvarezRusoSalaCOHGammaPXSec", pINFO)  << "*** Not in allowed kinematical space";
    return false;
  }

  // here we can add cuts on gamma energy to avoid peaks at maximum energy which are clearly unphysical  

  return true;
}
//____________________________________________________________________________
double 
AlvarezRusoSalaCOHGammaPXSec::NeutrinoHadronContraction( const Interaction & i, 
							 const utils::math::GTraceContraction & H ) const { 

  double k0 = i.InitState().ProbeE(kRfLab); // incident neutrino E
  TLorentzVector probe( 0., 0., k0, k0 ) ;
  TLorentzVector out_nu = i.Kine().FSLeptonP4();
  TLorentzVector t_q = probe - out_nu; // in/out neutrino E difference
  // FIXME not sure if q is exactly right
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() }; // Z boson momentum??

  std::complex<double> lh = -8*k0*((2*k0 - q[0] - q[3])*H(0,0,0,0) + q[1]*H(0,0,0,1) - std::complex<double>(0,1)*q[1]*H(0,0,0,2) +
				   (-2*k0 + q[0] + q[3])*H(0,0,0,3) + q[1]*H(0,1,0,0) + (-q[0] + q[3])*H(0,1,0,1) +
				   std::complex<double>(0,1)*(q[0] - q[3])*H(0,1,0,2) - q[1]*H(0,1,0,3) + std::complex<double>(0,1)*q[1]*H(0,2,0,0) -
				   std::complex<double>(0,1)*(q[0] - q[3])*H(0,2,0,1) + (-q[0] + q[3])*H(0,2,0,2) - std::complex<double>(0,1)*q[1]*H(0,2,0,3) +
				   (-2*k0 + q[0] + q[3])*H(0,3,0,0) - q[1]*H(0,3,0,1) + std::complex<double>(0,1)*q[1]*H(0,3,0,2) +
				   (2*k0 - q[0] - q[3])*H(0,3,0,3) + (-2*k0 + q[0] + q[3])*H(1,0,1,0) - q[1]*H(1,0,1,1) +
				   std::complex<double>(0,1)*q[1]*H(1,0,1,2) + (2*k0 - q[0] - q[3])*H(1,0,1,3) - q[1]*H(1,1,1,0) +
				   (q[0] - q[3])*H(1,1,1,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(1,1,1,2) + q[1]*H(1,1,1,3) -
				   std::complex<double>(0,1)*q[1]*H(1,2,1,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(1,2,1,1) + (q[0] - q[3])*H(1,2,1,2) +
				   std::complex<double>(0,1)*q[1]*H(1,2,1,3) + (2*k0 - q[0] - q[3])*H(1,3,1,0) + q[1]*H(1,3,1,1) -
				   std::complex<double>(0,1)*q[1]*H(1,3,1,2) + (-2*k0 + q[0] + q[3])*H(1,3,1,3) + (-2*k0 + q[0] + q[3])*H(2,0,2,0) -
				   q[1]*H(2,0,2,1) + std::complex<double>(0,1)*q[1]*H(2,0,2,2) + (2*k0 - q[0] - q[3])*H(2,0,2,3) - q[1]*H(2,1,2,0) +
				   (q[0] - q[3])*H(2,1,2,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(2,1,2,2) + q[1]*H(2,1,2,3) -
				   std::complex<double>(0,1)*q[1]*H(2,2,2,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(2,2,2,1) + (q[0] - q[3])*H(2,2,2,2) +
				   std::complex<double>(0,1)*q[1]*H(2,2,2,3) + (2*k0 - q[0] - q[3])*H(2,3,2,0) + q[1]*H(2,3,2,1) -
				   std::complex<double>(0,1)*q[1]*H(2,3,2,2) + (-2*k0 + q[0] + q[3])*H(2,3,2,3) + (-2*k0 + q[0] + q[3])*H(3,0,3,0) -
				   q[1]*H(3,0,3,1) + std::complex<double>(0,1)*q[1]*H(3,0,3,2) + (2*k0 - q[0] - q[3])*H(3,0,3,3) - q[1]*H(3,1,3,0) +
				   (q[0] - q[3])*H(3,1,3,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(3,1,3,2) + q[1]*H(3,1,3,3) -
				   std::complex<double>(0,1)*q[1]*H(3,2,3,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(3,2,3,1) + (q[0] - q[3])*H(3,2,3,2) +
				   std::complex<double>(0,1)*q[1]*H(3,2,3,3) + (2*k0 - q[0] - q[3])*H(3,3,3,0) + q[1]*H(3,3,3,1) -
				   std::complex<double>(0,1)*q[1]*H(3,3,3,2) + (-2*k0 + q[0] + q[3])*H(3,3,3,3));
  return lh.real();
}
//____________________________________________________________________________
double 
AlvarezRusoSalaCOHGammaPXSec::AntiNeutrinoHadronContraction( const Interaction & i, 
							     const utils::math::GTraceContraction & H ) const {

  double k0 = i.InitState().ProbeE(kRfLab); // incident neutrino E
  TLorentzVector probe( 0., 0., k0, k0 ) ;
  TLorentzVector out_nu = i.Kine().FSLeptonP4();
  TLorentzVector t_q = probe - out_nu; // in/out neutrino E difference
  // FIXME not sure if q is exactly right
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() }; // Z boson momentum??

  std::complex<double> lh = -8*k0*((2*k0 - q[0] - q[3])*H(0,0,0,0) + q[1]*H(0,0,0,1) +  std::complex<double>(0,1)*q[1]*H(0,0,0,2) +
				   (-2*k0 + q[0] + q[3])*H(0,0,0,3) + q[1]*H(0,1,0,0) + (-q[0] + q[3])*H(0,1,0,1) -
				   std::complex<double>(0,1)*(q[0] - q[3])*H(0,1,0,2) - q[1]*H(0,1,0,3) - std::complex<double>(0,1)*q[1]*H(0,2,0,0) +
				   std::complex<double>(0,1)*(q[0] - q[3])*H(0,2,0,1) + (-q[0] + q[3])*H(0,2,0,2) + std::complex<double>(0,1)*q[1]*H(0,2,0,3) +
				   (-2*k0 + q[0] + q[3])*H(0,3,0,0) - q[1]*H(0,3,0,1) - std::complex<double>(0,1)*q[1]*H(0,3,0,2) +
				   (2*k0 - q[0] - q[3])*H(0,3,0,3) + (-2*k0 + q[0] + q[3])*H(1,0,1,0) - q[1]*H(1,0,1,1) -
				   std::complex<double>(0,1)*q[1]*H(1,0,1,2) + (2*k0 - q[0] - q[3])*H(1,0,1,3) - q[1]*H(1,1,1,0) +
				   (q[0] - q[3])*H(1,1,1,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(1,1,1,2) + q[1]*H(1,1,1,3) +
				   std::complex<double>(0,1)*q[1]*H(1,2,1,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(1,2,1,1) + (q[0] - q[3])*H(1,2,1,2) -
				   std::complex<double>(0,1)*q[1]*H(1,2,1,3) + (2*k0 - q[0] - q[3])*H(1,3,1,0) + q[1]*H(1,3,1,1) +
				   std::complex<double>(0,1)*q[1]*H(1,3,1,2) + (-2*k0 + q[0] + q[3])*H(1,3,1,3) + (-2*k0 + q[0] + q[3])*H(2,0,2,0) -
				   q[1]*H(2,0,2,1) - std::complex<double>(0,1)*q[1]*H(2,0,2,2) + (2*k0 - q[0] - q[3])*H(2,0,2,3) - q[1]*H(2,1,2,0) +
				   (q[0] - q[3])*H(2,1,2,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(2,1,2,2) + q[1]*H(2,1,2,3) +
				   std::complex<double>(0,1)*q[1]*H(2,2,2,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(2,2,2,1) + (q[0] - q[3])*H(2,2,2,2) -
				   std::complex<double>(0,1)*q[1]*H(2,2,2,3) + (2*k0 - q[0] - q[3])*H(2,3,2,0) + q[1]*H(2,3,2,1) +
				   std::complex<double>(0,1)*q[1]*H(2,3,2,2) + (-2*k0 + q[0] + q[3])*H(2,3,2,3) + (-2*k0 + q[0] + q[3])*H(3,0,3,0) -
				   q[1]*H(3,0,3,1) - std::complex<double>(0,1)*q[1]*H(3,0,3,2) + (2*k0 - q[0] - q[3])*H(3,0,3,3) - q[1]*H(3,1,3,0) +
				   (q[0] - q[3])*H(3,1,3,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(3,1,3,2) + q[1]*H(3,1,3,3) +
				   std::complex<double>(0,1)*q[1]*H(3,2,3,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(3,2,3,1) + (q[0] - q[3])*H(3,2,3,2) -
				   std::complex<double>(0,1)*q[1]*H(3,2,3,3) + (2*k0 - q[0] - q[3])*H(3,3,3,0) + q[1]*H(3,3,3,1) +
				   std::complex<double>(0,1)*q[1]*H(3,3,3,2) + (-2*k0 + q[0] + q[3])*H(3,3,3,3));
  return lh.real();
}
//____________________________________________________________________________
void AlvarezRusoSalaCOHGammaPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSalaCOHGammaPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSalaCOHGammaPXSec::LoadConfig(void)
{

  fCurrents.clear() ;
  bool good_configuration = true ;

  // load all the resonance specific currents
  std::vector<RgKey> current_keys ;
  GetParamVectKeys( "COH-HadronicCurrents", current_keys ) ;

  for ( unsigned int i = 0 ; i < current_keys.size() ; ++i ) {

    const Algorithm * algo = SubAlg( current_keys[i] ) ;

    const COHHadronicCurrentI * current = dynamic_cast< const COHHadronicCurrentI * >( algo ) ;

    if ( ! current ) {
      good_configuration = false ;
      LOG("AlvarezRusoSalaCOHGammaPXSec", pERROR ) << "SubAlgo with key " << current_keys[i] << " not retrieved" ;
    }
    else {
      fCurrents.push_back( current ) ;
    }

  }


  //-- load the form factor
  fFormFactors =
    dynamic_cast<const COHFormFactorI *> (this->SubAlg("COH-FormFactor"));
  if (! fFormFactors ) {
    good_configuration = false ;
    LOG("AlvarezRusoSalaCOHGammaPXSec", pERROR ) << "Form factor not retrieved" ;
  }


  //-- load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  if ( ! fXSecIntegrator ) {
    good_configuration = false ;
    LOG("AlvarezRusoSalaCOHGammaPXSec", pERROR ) << "Integrator not retrieved" ;
  }

  if ( ! good_configuration ) {
    LOG("AlvarezRusoSalaCOHGammaPXSec", pFATAL ) << "Bad configuration: exiting" ;
    exit( 78 ) ;
  }
}
