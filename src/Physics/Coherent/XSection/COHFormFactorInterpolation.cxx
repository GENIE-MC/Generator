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
  COHFormFactorMap("genie::COHFormFactorInterpolation")
{

}
//____________________________________________________________________________
COHFormFactorInterpolation::COHFormFactorInterpolation(string config) :
  COHFormFactorMap("genie::COHFormFactorInterpolation", config)
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


}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & InterpolateProtons( int pdg ) const {

  return fInterProtons[pdg] = LinearInterpolation( pdg, pdg::IonPdgCodeToZ ) ;
}
//____________________________________________________________________________
const genie::FourierBesselFFCalculator & InterpolateNeutrons( int pdg ) const {

  return fInterNeutrons[pdg] =
    LinearInterpolation( pdg,
                        [](int pdg){ return pdg::IonPdgCodeToA(pdg)
                                            - pdg::IonPdgCodeToZ(pdg) } ) ;
}
//____________________________________________________________________________
genie::FourierBesselFFCalculator LinearInterpolation( int pdg,
                                                      const std::function<int(int)> & var ) {

    std::pair<int, int> neighbours = NearbyNuclei( pdg ) ;
    double r = RadiusInterpolation(pdg, neighbours) ;

    std::vector<double> first_v( Map()[neighbours.first] -> Coefficients() ) ;
    std::vector<double> second_v( Map()[neighbours.second] -> Coefficients() ) ;

    unsigned int n_coeffs = max( first_v.size(), second_v.size() ) ;

    // vector are extended with 0s if necessary
    first_v.resise( n_coeffs, 0. );
    second_v.resise( n_coeffs, 0. );

    vector<double> coeffs( 0., n_coeffs ) ;
    for ( unsinged int i = 0; i < n_coeffs; ++i ) {
      int var_first = var( neighbours.first) ;
      int var_second = var( neighbours.second ) ;

      if ( var_first == var_second ) {
        // in this case liner interpolation fails, we take the average
        coeffs[i] = 0.5 *( first_v[i] + second_v[i] );
      }
      else {
        int new_var = var(pdg) ;
        coeffs[i] = first_v[i] + (new_var - var_first)*(second_v[i] - first_v[i])/(var_second-var_first) ;
      }
    }

    return FourierBesselFFCalculator( coeffs, r ) ;
}

//____________________________________________________________________________
void COHFormFactorInterpolation::LoadConfig(void)
{

  COHFormFactorMap::LoadConfig() ;

  if ( ! good_configuration ) {
    LOG("COHFormFactorInterpolation", pFATAL ) << "Configuration not good, exiting" ;
    exit ( 78 ) ;
  }

}
//____________________________________________________________________________
