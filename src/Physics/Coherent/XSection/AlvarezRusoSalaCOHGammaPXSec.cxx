//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
 Author: Marco Roda
         University of Liverpool
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
  {

  }
//____________________________________________________________________________
AlvarezRusoSalaCOHGammaPXSec::AlvarezRusoSalaCOHGammaPXSec(string config) :
  XSecAlgorithmI("genie::AlvarezRusoSalaCOHGammaPXSec", config),
  fFormFactors( nullptr ),
  fXSecIntegrator( nullptr )
  {

  }
//____________________________________________________________________________
AlvarezRusoSalaCOHGammaPXSec::~AlvarezRusoSalaCOHGammaPXSec()
{

}
//____________________________________________________________________________
double AlvarezRusoSalaCOHGammaPXSec::XSec( const Interaction * interaction,
                                           KinePhaseSpace_t kps) const {

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  int pdg = init_state.Tgt().Pdg() ;
  double E_nu = init_state.ProbeE(kRfLab); // neutrino energy

  const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  const TLorentzVector p4_gam = kinematics.HadSystP4();
  double E_lep = p4_lep.E();

  // get the traces for each currents and sum them
  double xsec = 0.0;


  return (xsec);
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
