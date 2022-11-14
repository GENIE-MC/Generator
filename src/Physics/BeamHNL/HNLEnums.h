//----------------------------------------------------------------------------
/*!

  All enums used in the HNL simulation (should) live here.

\namespace  genie::HNL::HNLenums

\brief      Typedef enums

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 10th, 2022

\cpright    Copyright (c) 2003-2022, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//----------------------------------------------------------------------------

#ifndef _HNL_ENUMS_H_
#define _HNL_ENUMS_H_

// -- C++ includes

#include <map>

// -- GENIE includes

#include "Framework/Conventions/Units.h"

namespace genie {
namespace HNL {

    namespace HNLenums {

	/// Coupling configuration indices
	typedef enum t_coupIdx {

	    kInitIdx     = -1,
	    kEqualIdx    = 0,
	    kMuonIdx     = 1,
	    kElectronIdx = 2,
	    kOtherIdx    = 3

	} coupIdx_t;

	typedef enum t_parent {

	    kNoPar = -1,
	    kAll   = 0,
	    kPion  = 1,
	    kKaon  = 2,
	    kMuon  = 3,
	    kNeuk  = 4
	    
	} parent_t;

	typedef enum t_nutype {

	  kNone = -1,
	  kNumu = 1,
	  kNumubar = 2,
	  kNue = 3,
	  kNuebar = 4

	} nutype_t;

    } // namespace HNLenums

} // namespace HNL
} // namespace genie

#endif // #ifndef _HNL_ENUMS_H_
