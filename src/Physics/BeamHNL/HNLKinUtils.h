//----------------------------------------------------------------------------
/*!

\namespace  genie::utils::hnl

\brief      Useful kinematic functions

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 11th, 2022

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _HNL_KINUTILS_H_
#define _HNL_KINUTILS_H_

// -- C++ includes
#include <cmath>

// -- GENIE includes
#include "Framework/Messenger/Messenger.h"

namespace genie{
namespace utils {

    namespace hnl {

	inline double MassX( double m1, double m2 ) {
	    if( m2 <= 0. || m1 < 0.) { LOG( "HNL", pERROR ) << "Illegal masses m1 = " << m1 << ", m2 = " << m2; exit( 3 ); }
	    return m1 / m2;
	}
	
	inline double Kallen( double x, double y, double z ) {
	    return x*x + y*y + z*z - 2. * ( x*y + y*z + z*x );
	}

	inline double SymmDiff( double x, double y ) {
	  return x + y - ( x-y ) * ( x-y );
	}

	inline double RhoFunc( double x, double y ) {
	  return SymmDiff( x, y ) * std::sqrt( Kallen( 1, x, y ) );
	}
	
    } // namespace hnl

} //namespace utils
} //namespace genie

#endif // #ifndef _HNL_KINUTILS_H_
