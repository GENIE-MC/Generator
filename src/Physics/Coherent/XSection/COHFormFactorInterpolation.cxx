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

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/COHFormFactorInterpolation.h"


using namespace genie;


COHFormFactorInterpolation::COHFormFactorInterpolation() :
  COHFormFactorMap("genie::COHFormFactorInterpolation"),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHFormFactorInterpolation::COHFormFactorInterpolation(string config) :
  COHFormFactorMap("genie::COHFormFactorInterpolation", config),
  fAllowExtrapolation(false)
{

}
//____________________________________________________________________________
COHFormFactorInterpolation::~COHFormFactorInterpolation()
{

}
//____________________________________________________________________________
double COHFormFactorInterpolation::ProtonFF( double Q, int pdg ) const {

  if ( ! HasNucleus(pdg) ) return 0. ;

  if ( COHFormFactorMap::HasNucleus(pdg) ) return COHFormFactorMap::ProtonFF( Q, pdg ) ;

  const std::map<int,genie::FourierBesselFFCalculator>::const_iterator it =
    fInterProtons.find( pdg ) ;

  if ( it != fInterProtons.end() )  return it.second.FormFactor(Q) ;

  return InterpolateProtons( pdg ).FormFactor(Q) ;

}
//____________________________________________________________________________
double COHFormFactorInterpolation::NeutronFF( double Q, int pdg ) const {

  if ( ! HasNucleus(pdg) ) return 0. ;

  if ( COHFormFactorMap::HasNucleus(pdg) ) return COHFormFactorMap::NeutronFF( Q, pdg ) ;

  const std::map<int,genie::FourierBesselFFCalculator>::const_iterator it =
    fInterNeutrons.find( pdg ) ;

  if ( it != fInterNeutrons.end() )  return it.second.FormFactor(Q) ;

  return InterpolateNeutrons( pdg ).FormFactor(Q) ;

}
//____________________________________________________________________________
bool COHFormFactorInterpolation::HasNucleus( int pdg ) const {

  if ( Map().size() == 0 ) return false ;

   if ( fAllowExtrapolation ) {
     if ( pdg >= Map().begin()->first ) return true ;
   }
   else {
     if ( pdg < Map().begin()->first ) return false ;
     else if ( pdg > Map().rbegin()->first ) return false ;
     else return true ;
   }

   return true ;
}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & COHFormFactorInterpolation::InterpolateProtons( int pdg ) const {

  return fInterProtons[pdg] = LinearInterpolation( pdg, pdg::IonPdgCodeToZ ) ;
}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & COHFormFactorInterpolation::InterpolateNeutrons( int pdg ) const {

  return fInterNeutrons[pdg] =
    LinearInterpolation( pdg,
                        [](int pdg){ return pdg::IonPdgCodeToA(pdg)
                                            - pdg::IonPdgCodeToZ(pdg) } ) ;
}
//____________________________________________________________________________
pair<int, int> COHFormFactorInterpolation::NearbyNuclei( pdg ) const {

  auto last = Map().rbegin() ;

  if ( last.first < pdg ) return make_pair( (last+1).first, last.first ) ;

  auto first = 1+last ;

  while (first.first > pdg ) {
    first++
  }

  return make_pair( first.first, (first+1).first ) ;
}
//____________________________________________________________________________
double COHFormFactorInterpolation::RadiusInterpolation( int pdg,
                                                        const pair<int, int> & neighbours ) const {

      // we interpolate the radius according to
      // r = r0 + c pow( A, 1/3 )
      // which is a linear interpolation as a function of  pow( A, 1/3 )

      int min_A = pdg::IonPdgCodeToA(neighbours.first) ;
      int max_A = pdg::IonPdgCodeToA(neighbours.second) ;

      double r = 0. ;
      if ( min_A == max_A ) {
        // take the average
        r = 0.5( Map()[neighbours.first].Calculator().Radius()
                 + Map()[neighbours.second].Calculator().Radius() ) ;
      }
      else {
        double min_root = pow( min_A, 1./3 ) ;
        double max_root = pow( max_A, 1./3 ) ;
        double root = pow( pdg::IonPdgCodeToA(pdg), 1./3 ) ;

        double min_rad =   Map()[neighbours.first].Calculator().Radius() ;
        double max_rad =   Map()[neighbours.second].Calculator().Radius() ;

        r = min_rad
            + (max_rad - min_rad)*(root - min_root)/(max_root - min_root) ;

      }

      return r ;
}

//____________________________________________________________________________

genie::FourierBesselFFCalculator COHFormFactorInterpolation::LinearInterpolation( int pdg,
                                                                                  const std::function<int(int)> & var ) {

    std::pair<int, int> neighbours = NearbyNuclei( pdg ) ;

    double r = RadiusInterpolation(pdg, neighbours) ;

    std::vector<double> first_v( Map()[neighbours.first] -> Coefficients() ) ;
    std::vector<double> second_v( Map()[neighbours.second] -> Coefficients() ) ;

    unsigned int n_coeffs = max( first_v.size(), second_v.size() ) ;

    // vector are extended with 0s if necessary
    first_v.resize( n_coeffs, 0. );
    second_v.resize( n_coeffs, 0. );

    vector<double> coeffs( 0., n_coeffs ) ;

    int var_first = var( neighbours.first) ;
    int var_second = var( neighbours.second ) ;

    if ( var_first == var_second ) {
      // in this case liner interpolation fails, we take the average
      for ( unsinged int i = 0; i < n_coeffs; ++i ) {
        coeffs[i] = 0.5 *( first_v[i] + second_v[i] );
      }
    }
    else {

      int new_var = var(pdg) ;
      double scale = (new_var - var_first)/ (double) (var_second-var_first) ;

      for ( unsinged int i = 0; i < n_coeffs; ++i ) {
        coeffs[i] = first_v[i] + scale*(second_v[i] - first_v[i])/ ;
      }
    }

    return FourierBesselFFCalculator( coeffs, r ) ;
}

//____________________________________________________________________________
void COHFormFactorInterpolation::LoadConfig(void)
{

  COHFormFactorMap::LoadConfig() ;

  bool good_configuration = true ;
  if ( Map().size() < 2 ) {
    LOG("COHFormFactorInterpolation", pERROR ) << "Not enough FormFactors inside the Map" ;
    good_configuration = false ;
  }

  GetParamDef( "AllowExtrapolation", fAllowExtrapolation, false ) ;

  if ( ! good_configuration ) {
    LOG("COHFormFactorInterpolation", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
