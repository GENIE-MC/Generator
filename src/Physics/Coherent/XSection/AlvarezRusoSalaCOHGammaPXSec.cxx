//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
<<<<<<< HEAD

 Author: Marco Roda
         University of Liverpool
         July 2020

 For the class documentation see the corresponding header file.

=======
 Author: Marco Roda
         University of Liverpool
         July 2020
 For the class documentation see the corresponding header file.
>>>>>>> eace98236f5cbc8f7823245bfd3ec71466eece6b
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
                                           KinePhaseSpace_t kps) const {

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  // int pdg = init_state.Tgt().Pdg() ;
  double E_nu = init_state.ProbeE(kRfLab); // neutrino energy

  // const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  const TLorentzVector p4_gam = kinematics.HadSystP4();
  // double E_lep = p4_lep.E();
  double E_gamma = p4_gam.E();

  // get the traces for each currents and sum them
  utils::math::GTrace R = HadronicCurrent(interaction);

  int nu = init_state.ProbePdg();
  double contraction = 0.;

  if (pdg::IsNeutrino(nu)) {
	contraction = NeutrinoHadronContraction(interaction, R);
  }
  else if (pdg::IsAntiNeutrino(nu)) {
    contraction = AntiNeutrinoHadronContraction(interaction, R);
  }

  double pre_factor = 4.*constants::kPi*constants::kAem * constants::kGF2 / 16*pow(2.*constants::kPi, 5);

  return pre_factor*contraction*( E_gamma*(E_nu - E_gamma) / E_nu );
}
//____________________________________________________________________________
utils::math::GTrace AlvarezRusoSalaCOHGammaPXSec::HadronicCurrent( const Interaction * interaction ) const {
  // FIXME
  const COHFormFactorI * ff;

  // sum R's from all currents here
  utils::math::GTrace R_total;
  for ( unsigned int i = 0; i < fCurrents.size(); i++ ) {
    R_total += fCurrents[i] -> R(interaction, ff);
  }

  return R_total;
}
//____________________________________________________________________________
std::complex<double> AlvarezRusoSalaCOHGammaPXSec::H( const utils::math::GTrace & R,  unsigned short i,
                                      unsigned short j, unsigned short k, unsigned short l ) const
{
	return R[i][j] * std::conj(R[k][l]);
>>>>>>> eace98236f5cbc8f7823245bfd3ec71466eece6b
}
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::Integral(const Interaction * interaction) const
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
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::AntiNeutrinoHadronContraction( const Interaction * i, const utils::math::GTrace & R ) const {

	TLorentzVector *probe = i -> InitState().GetProbeP4( kRfLab );
	TLorentzVector out_nu = i -> Kine().FSLeptonP4();
	TLorentzVector t_q = *probe - out_nu;
	// FIXME not sure if q is exactly right
	std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() };
	double k0 = *probe->E();
	delete probe;


	std::complex<double> lh = -8*k0*((2*k0 - q[0] - q[3])*H(R,0,0,0,0) + q[1]*H(R,0,0,0,1) - std::complex<double>(0,1)*q[1]*H(R,0,0,0,2) +
	                       (-2*k0 + q[0] + q[3])*H(R,0,0,0,3) + q[1]*H(R,0,1,0,0) + (-q[0] + q[3])*H(R,0,1,0,1) +
	                       std::complex<double>(0,1)*(q[0] - q[3])*H(R,0,1,0,2) - q[1]*H(R,0,1,0,3) + std::complex<double>(0,1)*q[1]*H(R,0,2,0,0) -
	                       std::complex<double>(0,1)*(q[0] - q[3])*H(R,0,2,0,1) + (-q[0] + q[3])*H(R,0,2,0,2) - std::complex<double>(0,1)*q[1]*H(R,0,2,0,3) +
	                       (-2*k0 + q[0] + q[3])*H(R,0,3,0,0) - q[1]*H(R,0,3,0,1) + std::complex<double>(0,1)*q[1]*H(R,0,3,0,2) +
	                       (2*k0 - q[0] - q[3])*H(R,0,3,0,3) + (-2*k0 + q[0] + q[3])*H(R,1,0,1,0) - q[1]*H(R,1,0,1,1) +
	                       std::complex<double>(0,1)*q[1]*H(R,1,0,1,2) + (2*k0 - q[0] - q[3])*H(R,1,0,1,3) - q[1]*H(R,1,1,1,0) +
	                       (q[0] - q[3])*H(R,1,1,1,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,1,1,1,2) + q[1]*H(R,1,1,1,3) -
	                       std::complex<double>(0,1)*q[1]*H(R,1,2,1,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,1,2,1,1) + (q[0] - q[3])*H(R,1,2,1,2) +
	                       std::complex<double>(0,1)*q[1]*H(R,1,2,1,3) + (2*k0 - q[0] - q[3])*H(R,1,3,1,0) + q[1]*H(R,1,3,1,1) -
	                       std::complex<double>(0,1)*q[1]*H(R,1,3,1,2) + (-2*k0 + q[0] + q[3])*H(R,1,3,1,3) + (-2*k0 + q[0] + q[3])*H(R,2,0,2,0) -
	                       q[1]*H(R,2,0,2,1) + std::complex<double>(0,1)*q[1]*H(R,2,0,2,2) + (2*k0 - q[0] - q[3])*H(R,2,0,2,3) - q[1]*H(R,2,1,2,0) +
	                       (q[0] - q[3])*H(R,2,1,2,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,2,1,2,2) + q[1]*H(R,2,1,2,3) -
	                       std::complex<double>(0,1)*q[1]*H(R,2,2,2,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,2,2,2,1) + (q[0] - q[3])*H(R,2,2,2,2) +
	                       std::complex<double>(0,1)*q[1]*H(R,2,2,2,3) + (2*k0 - q[0] - q[3])*H(R,2,3,2,0) + q[1]*H(R,2,3,2,1) -
	                       std::complex<double>(0,1)*q[1]*H(R,2,3,2,2) + (-2*k0 + q[0] + q[3])*H(R,2,3,2,3) + (-2*k0 + q[0] + q[3])*H(R,3,0,3,0) -
	                       q[1]*H(R,3,0,3,1) + std::complex<double>(0,1)*q[1]*H(R,3,0,3,2) + (2*k0 - q[0] - q[3])*H(R,3,0,3,3) - q[1]*H(R,3,1,3,0) +
	                       (q[0] - q[3])*H(R,3,1,3,1) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,3,1,3,2) + q[1]*H(R,3,1,3,3) -
	                       std::complex<double>(0,1)*q[1]*H(R,3,2,3,0) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,3,2,3,1) + (q[0] - q[3])*H(R,3,2,3,2) +
	                       std::complex<double>(0,1)*q[1]*H(R,3,2,3,3) + (2*k0 - q[0] - q[3])*H(R,3,3,3,0) + q[1]*H(R,3,3,3,1) -
	                       std::complex<double>(0,1)*q[1]*H(R,3,3,3,2) + (-2*k0 + q[0] + q[3])*H(R,3,3,3,3));
  return lh.real();
}
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::NeutrinoHadronContraction( const Interaction * i, const utils::math::GTrace & R ) const {

	TLorentzVector *probe = i -> InitState().GetProbeP4( kRfLab );
	TLorentzVector out_nu = i -> Kine().FSLeptonP4();
	TLorentzVector t_q = *probe - out_nu;
	// FIXME not sure if q is exactly right
	std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() };
	double k0 = *probe->E();
	delete probe;

	std::complex<double> lh = -8*k0*((2*k0 - q[0] - q[3])*H(R,0,0,0,0) + q[1]*H(R,0,0,0,1) +  std::complex<double>(0,1)*q[1]*H(R,0,0,0,2) +
							 (-2*k0 + q[0] + q[3])*H(R,0,0,0,3) + q[1]*H(R,0,1,0,0) + (-q[0] + q[3])*H(R,0,1,0,1) -
							 std::complex<double>(0,1)*(q[0] - q[3])*H(R,0,1,0,2) - q[1]*H(R,0,1,0,3) - std::complex<double>(0,1)*q[1]*H(R,0,2,0,0) +
							 std::complex<double>(0,1)*(q[0] - q[3])*H(R,0,2,0,1) + (-q[0] + q[3])*H(R,0,2,0,2) + std::complex<double>(0,1)*q[1]*H(R,0,2,0,3) +
							 (-2*k0 + q[0] + q[3])*H(R,0,3,0,0) - q[1]*H(R,0,3,0,1) - std::complex<double>(0,1)*q[1]*H(R,0,3,0,2) +
							 (2*k0 - q[0] - q[3])*H(R,0,3,0,3) + (-2*k0 + q[0] + q[3])*H(R,1,0,1,0) - q[1]*H(R,1,0,1,1) -
							 std::complex<double>(0,1)*q[1]*H(R,1,0,1,2) + (2*k0 - q[0] - q[3])*H(R,1,0,1,3) - q[1]*H(R,1,1,1,0) +
							 (q[0] - q[3])*H(R,1,1,1,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,1,1,1,2) + q[1]*H(R,1,1,1,3) +
							 std::complex<double>(0,1)*q[1]*H(R,1,2,1,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,1,2,1,1) + (q[0] - q[3])*H(R,1,2,1,2) -
							 std::complex<double>(0,1)*q[1]*H(R,1,2,1,3) + (2*k0 - q[0] - q[3])*H(R,1,3,1,0) + q[1]*H(R,1,3,1,1) +
							 std::complex<double>(0,1)*q[1]*H(R,1,3,1,2) + (-2*k0 + q[0] + q[3])*H(R,1,3,1,3) + (-2*k0 + q[0] + q[3])*H(R,2,0,2,0) -
							 q[1]*H(R,2,0,2,1) - std::complex<double>(0,1)*q[1]*H(R,2,0,2,2) + (2*k0 - q[0] - q[3])*H(R,2,0,2,3) - q[1]*H(R,2,1,2,0) +
							 (q[0] - q[3])*H(R,2,1,2,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,2,1,2,2) + q[1]*H(R,2,1,2,3) +
							 std::complex<double>(0,1)*q[1]*H(R,2,2,2,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,2,2,2,1) + (q[0] - q[3])*H(R,2,2,2,2) -
							 std::complex<double>(0,1)*q[1]*H(R,2,2,2,3) + (2*k0 - q[0] - q[3])*H(R,2,3,2,0) + q[1]*H(R,2,3,2,1) +
							 std::complex<double>(0,1)*q[1]*H(R,2,3,2,2) + (-2*k0 + q[0] + q[3])*H(R,2,3,2,3) + (-2*k0 + q[0] + q[3])*H(R,3,0,3,0) -
							 q[1]*H(R,3,0,3,1) - std::complex<double>(0,1)*q[1]*H(R,3,0,3,2) + (2*k0 - q[0] - q[3])*H(R,3,0,3,3) - q[1]*H(R,3,1,3,0) +
							 (q[0] - q[3])*H(R,3,1,3,1) + std::complex<double>(0,1)*(q[0] - q[3])*H(R,3,1,3,2) + q[1]*H(R,3,1,3,3) +
							 std::complex<double>(0,1)*q[1]*H(R,3,2,3,0) - std::complex<double>(0,1)*(q[0] - q[3])*H(R,3,2,3,1) + (q[0] - q[3])*H(R,3,2,3,2) -
							 std::complex<double>(0,1)*q[1]*H(R,3,2,3,3) + (2*k0 - q[0] - q[3])*H(R,3,3,3,0) + q[1]*H(R,3,3,3,1) +
							 std::complex<double>(0,1)*q[1]*H(R,3,3,3,2) + (-2*k0 + q[0] + q[3])*H(R,3,3,3,3));


  return lh.real();

}
