//____________________________________________________________________________
/*
  Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <vector>
#include <utility>
#include <algorithm>

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/COHProtonFormFactorInterpolation.h"


using namespace genie;


COHProtonFormFactorInterpolation::COHProtonFormFactorInterpolation() :
  COHFormFactorI("genie::COHProtonFormFactorInterpolation", "Default" ),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHProtonFormFactorInterpolation::COHProtonFormFactorInterpolation(string config) :
  COHFormFactorI("genie::COHProtonFormFactorInterpolation", config),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHProtonFormFactorInterpolation::~COHProtonFormFactorInterpolation()
{

}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::ProtonFF( double Q, int pdg ) const {

  if ( ! HasNucleus(pdg) ) return 0. ;

  if ( DeVriesFormFactorMap::HasNucleus(pdg) ) return DeVriesFormFactorMap::ProtonFF( Q, pdg ) ;

  const std::map<int,genie::FourierBesselFFCalculator>::const_iterator it =
    fInterProtons.find( pdg ) ;

  if ( it != fInterProtons.end() )  return it -> second.FormFactor(Q) ;

  return InterpolateProtons( pdg ).FormFactor(Q) ;

}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::NeutronFF( double Q, int pdg ) const {

  int z = pdg::IonPdgCodeToZ( pdg ) ; 
  
  double scale = ( pdg::IonPdgCodeToA( pdg ) - z ) / (double) z ; 

  return scale * ProtonFF( Q, pdg ) ; 
  
}
//____________________________________________________________________________
bool COHProtonFormFactorInterpolation::HasNucleus( int pdg ) const {

  if ( Map().size() == 0 ) return false ;

   if ( fAllowExtrapolation ) {
     if ( pdg >= Map().begin()->first ) return true ;
     else return false;
   }
   else {
     if ( pdg < Map().begin()->first ) return false ;
     else if ( pdg > Map().rbegin()->first ) return false ;
     else return true ;
   }

   return true ;
}
//____________________________________________________________________________
genie::Range1D_t COHProtonFormFactorInterpolation::QRange( int pdg ) const { 

  if ( ! HasNucleus(pdg) ) return COHFormFactorI::QRange( pdg ) ;

  if ( DeVriesFormFactorMap::HasNucleus(pdg) ) return DeVriesFormFactorMap::QRange( pdg ) ;

  std::map<int,genie::FourierBesselFFCalculator>::const_iterator proton_it =
    fInterProtons.find( pdg ) ;

  if ( proton_it == fInterProtons.end() ) {
    InterpolateProtons( pdg ) ;
    proton_it = fInterProtons.find( pdg ) ;
  }

  std::map<int,genie::FourierBesselFFCalculator>::const_iterator neutron_it =
    fInterNeutrons.find( pdg ) ;
  
  if ( neutron_it == fInterNeutrons.end() ) {
    InterpolateNeutrons( pdg ) ;
    neutron_it = fInterNeutrons.find( pdg ) ;
  }

  return Range1D_t( TMath::Min( proton_it -> second.QMin(), neutron_it -> second.QMin() ), 
		    TMath::Max( proton_it -> second.QMax(), neutron_it -> second.QMax() ) ) ;
		      
}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & COHProtonFormFactorInterpolation::InterpolateProtons( int pdg ) const {

  auto result = fInterProtons.insert( std::make_pair( pdg, LinearInterpolation( pdg, pdg::IonPdgCodeToZ ) ) ) ;
  return result.first -> second ;

}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & COHProtonFormFactorInterpolation::InterpolateNeutrons( int pdg ) const {

  auto result = fInterNeutrons.insert( std::make_pair( pdg, 
						       LinearInterpolation( pdg, 
									    [](int _pdg){ return pdg::IonPdgCodeToA(_pdg) 
                                                                                           - pdg::IonPdgCodeToZ(_pdg); }) ) ) ;
  return result.first -> second ;
}
//____________________________________________________________________________
pair<int, int> COHProtonFormFactorInterpolation::NearbyNuclei( int pdg ) const {

  auto rit = Map().rbegin() ;

  if ( rit -> first < pdg ) return std::make_pair( (++Map().rbegin()) -> first, 
						   rit -> first ) ;

  do {
   rit++ ;
  } while ( rit -> first > pdg && rit != Map().rend() ) ;
 
  auto sec = rit ;
  sec-- ;
  return std::make_pair( rit -> first, sec -> first ) ;
}
//____________________________________________________________________________
double COHProtonFormFactorInterpolation::RadiusInterpolation( int pdg,
                                                        const pair<int, int> & neighbours ) const {

      // we interpolate the radius according to
      // r = r0 + c pow( A, 1/3 )
      // which is a linear interpolation as a function of  pow( A, 1/3 )

      int min_A = pdg::IonPdgCodeToA(neighbours.first) ;
      int max_A = pdg::IonPdgCodeToA(neighbours.second) ;

      double r = 0. ;
      if ( min_A == max_A ) {
        // take the average
        r = 0.5 * ( Map().at(neighbours.first) -> Calculator().Radius()
                    + Map().at(neighbours.second) -> Calculator().Radius() ) ;
      }
      else {
        double min_root = pow( min_A, 1./3 ) ;
        double max_root = pow( max_A, 1./3 ) ;
        double root = pow( pdg::IonPdgCodeToA(pdg), 1./3 ) ;

        double min_rad =   Map().at( neighbours.first ) -> Calculator().Radius() ;
        double max_rad =   Map().at( neighbours.second ) -> Calculator().Radius() ;

        r = min_rad
            + (max_rad - min_rad)*(root - min_root)/(max_root - min_root) ;

      }

      return r ;
}

//____________________________________________________________________________

genie::FourierBesselFFCalculator COHProtonFormFactorInterpolation::LinearInterpolation( int pdg,
                                                                                  const std::function<int(int)> & var ) const {

    std::pair<int, int> neighbours = NearbyNuclei( pdg ) ;

    double r = RadiusInterpolation(pdg, neighbours) ;

    const genie::FourierBesselFFCalculator & first_calc = Map().at( neighbours.first) -> Calculator() ;
    const genie::FourierBesselFFCalculator & second_calc = Map().at( neighbours.second ) -> Calculator() ;
    
    std::vector<double> first_v( first_calc.Coefficients() ) ;
    std::vector<double> second_v( second_calc.Coefficients() ) ;

    unsigned int n_coeffs = std::max( first_v.size(), second_v.size() ) ;

    // vector are extended with 0s if necessary
    first_v.resize( n_coeffs, 0. );
    second_v.resize( n_coeffs, 0. );

    vector<double> coeffs(n_coeffs, 0. ) ;

    int var_first = var( neighbours.first ) ;
    int var_second = var( neighbours.second ) ;

    std::pair<double, double> range_first ( first_calc.QMin(), first_calc.QMax() ) ;
    std::pair<double, double> range_second ( first_calc.QMin(), first_calc.QMax() ) ;
    std::pair<double, double> new_range ;    

    if ( var_first == var_second ) {
      // in this case liner interpolation fails, we take the average
      for ( unsigned int i = 0; i < n_coeffs; ++i ) {
        coeffs[i] = 0.5 *( first_v[i] + second_v[i] );
      }
      new_range = std::make_pair( 0.5*( range_first.first +  range_second.first ), 
				  0.5*( range_first.second + range_second.second ) ) ;
    }
    else {

      int new_var = var(pdg) ;
      double scale = (new_var - var_first)/ (double) (var_second-var_first) ;

      for ( unsigned int i = 0; i < n_coeffs; ++i ) {
        coeffs[i] = first_v[i] + scale*(second_v[i] - first_v[i]) ;
      }

      new_range = std::make_pair( range_first.first  + scale * ( range_second.first  - range_first.first ), 
				  range_first.second + scale * ( range_second.second - range_first.second ) ) ; 

    }

    return FourierBesselFFCalculator( coeffs, r , new_range.first, new_range.second ) ;
}

//____________________________________________________________________________
void COHProtonFormFactorInterpolation::LoadConfig(void)
{

  DeVriesFormFactorMap::LoadConfig() ;

  bool good_configuration = true ;
  if ( Map().size() < 2 ) {
    LOG("COHProtonFormFactorInterpolation", pERROR ) << "Not enough FormFactors inside the Map" ;
    good_configuration = false ;
  }

  GetParamDef( "AllowExtrapolation", fAllowExtrapolation, false ) ;

  if ( ! good_configuration ) {
    LOG("COHProtonFormFactorInterpolation", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
