//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cmath>
#include <limits>
#include <algorithm>


#include "Physics/Coherent/XSection/COHGammaIntegrationLimits.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Controls.h" 
#include "Framework/Conventions/Constants.h" 


using namespace genie;


//____________________________________________________________________________
COHGammaIntegrationLimits::COHGammaIntegrationLimits() :
  Algorithm("genie::COHGammaIntegrationLimits"), 
  fFF( nullptr ), 
  fMaxEg( std::numeric_limits<double>::infinity() ) 
{ ; }
//____________________________________________________________________________
COHGammaIntegrationLimits::COHGammaIntegrationLimits(string config) :
  Algorithm("genie::COHGammaIntegrationLimits", config),
  fFF( nullptr ),
  fMaxEg( std::numeric_limits<double>::infinity() ) 
 { ; }
//____________________________________________________________________________
COHGammaIntegrationLimits::~COHGammaIntegrationLimits()
{

}
//____________________________________________________________________________
Range1D_t COHGammaIntegrationLimits::EGamma( const Interaction & in ) const {

  return Range1D_t( controls::kASmallNum, 
		    std::min( in.InitState().ProbeE( kRfLab ) - controls::kASmallNum, 
			      fMaxEg ) 
		    ) ; 
}
//____________________________________________________________________________
Range1D_t COHGammaIntegrationLimits::ThetaGamma( const Interaction & ) const {

  return Range1D_t( controls::kASmallNum, 
		    constants::kPi - controls::kASmallNum ) ; 
}
//____________________________________________________________________________
Range1D_t COHGammaIntegrationLimits::PhiGamma( const Interaction & ) const {

  return Range1D_t( 0., 
		    2*constants::kPi - controls::kASmallNum ) ; 
}
//____________________________________________________________________________
Range1D_t COHGammaIntegrationLimits::ThetaLepton( const Interaction & in ) const {

  double e_nu = in.InitState().ProbeE( kRfLab ) ; 
  double max_t = t( in ).max ;
  Range1D_t e_gamma_range = EGamma( in ) ; 
  double target_mass = in.InitState().Tgt().Mass() ;

  double min_e_l = e_nu - e_gamma_range.max - 0.5 * max_t / target_mass ;

  double min = controls::kASmallNum ; 
  double max = constants::kPi - controls::kASmallNum ; 

  if ( min_e_l > 0. ) {

    double sqrt_s = sqrt( target_mass * target_mass + 2. * e_nu * target_mass ) ; 
    
    double min_cos_theta_limit = 1. - 2. * target_mass * ( sqrt_s - target_mass - e_gamma_range.min ) / ( sqrt_s * min_e_l ) ; 
    
    double min_cos_theta = std::max( -1., min_cos_theta_limit ) ;

    max = acos( min_cos_theta ) ;
  }

  return Range1D_t( min, max ) ; 
}
//____________________________________________________________________________
Range1D_t COHGammaIntegrationLimits::t( const Interaction & i ) const {

  Range1D_t Q_range = fFF -> QRange( i.InitState().Tgt().Pdg() ) ;
  return Range1D_t( Q_range.min*Q_range.min , 
		    Q_range.max*Q_range.max ) ;

}
//____________________________________________________________________________
void COHGammaIntegrationLimits::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHGammaIntegrationLimits::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHGammaIntegrationLimits::LoadConfig(void)
{

  bool good_configuration = true ;

  //-- load the form factor                                                                                          
  fFF = dynamic_cast<const COHFormFactorI *> (this->SubAlg("COH-FormFactor"));
  if (! fFF ) {
    good_configuration = false ;
    LOG("COHGammaIntegrationLimits", pERROR ) << "Form factor not retrieved" ;
  }

  GetParam( "MaxGammaEnergy", fMaxEg ) ;

  if ( ! good_configuration ) {
    LOG("COHGammaIntegrationLimits", pFATAL ) << "Bad configuration: exiting" ;
    exit( 78 ) ;
  }


}
//____________________________________________________________________________
