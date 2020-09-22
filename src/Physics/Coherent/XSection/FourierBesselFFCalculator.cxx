//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda
         University of Liverpool

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cmath>

#include "Physics/Coherent/XSection/FourierBesselFFCalculator.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

using namespace genie;

//____________________________________________________________________________
double FourierBesselFFCalculator::FormFactor( double Q ) const {

  // Only keep Q values within valid range for parameters
  if ( Q > fQmax ) return 0. ;
  if ( Q < fQmin ) return 0. ;

  double qr = Q * fRadius ;

  double aux_sum = 0.0, nu ;

  for ( int i = fFBCs.size() - 1 ; i >= 0 ; --i ) {
     nu = i + 1. ;
     double pi_x_i = constants::kPi*nu ;
     aux_sum += pow( -1.0, i )*fFBCs[i]/( ( pi_x_i + qr )*( pi_x_i - qr ) ) ;
 }

 return 4.*constants::kPi*pow( fRadius/units::fm, 3)*aux_sum*(sin(qr)/(qr) ) ;

}
//____________________________________________________________________________
