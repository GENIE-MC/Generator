//----------------------------------------------------------------------------
/*!

\namespace  genie::utils::nhl

\brief      Useful kinematic functions

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 11th, 2022

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------

#ifndef _NHL_KINUTILS_H_
#define _NHL_KINUTILS_H_

// -- GENIE includes
#include "Framework/Messenger/Messenger.h"

namespace genie{
namespace utils {

    namespace nhl {

	inline double MassX( double m1, double m2 ) {
	    if( m2 <= 0. || m1 < 0.) { LOG( "SimpleHNL", pERROR ) << "BRFunctions::MassX:: Illegal masses m1 = " << m1 << ", m2 = " << m2; exit( 3 ); }
	    return m1 / m2;
	}
	
	inline double Kallen( double x, double y, double z ) {
	    return x*x + y*y + z*z - 2. * ( x*y + y*z + z*x );
	}
	
    } // namespace nhl

} //namespace utils
} //namespace genie

#endif // #ifndef _NHL_KINUTILS_H_
