//____________________________________________________________________________
/*!

\class    genie::hnl::HNLDecayMode

\brief    Enumeration of HNL decay modes.

\author   John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford

	  Costas Andreopoulos <c.andreopoulos \at cern.ch>
	  University of Liverpool

\created  November 10, 2011

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _HNL_DECAY_MODE_H_
#define _HNL_DECAY_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace hnl {

 typedef enum EHNLDecayMode {

   kHNLDcyNull     = -1, // dummy
   kHNLDcyNuNuNu   = 0,	 // N --> 3 nus. Summed over all flavours 
   kHNLDcyNuEE     = 1,	 // N --> nu_{a}      e^{\mp}   e^{\pm}. W and Z interfere
   kHNLDcyNuMuE    = 2,	 // N --> nu_{\mu/e}  e^{\mp} \mu^{\pm}. Only W. Summed over nue and numu
   kHNLDcyPi0Nu    = 3,	 // N --> \pi^{0}   \nu (any kind)	  
   kHNLDcyPiE      = 4,  // N --> \pi^{\pm}   e^{\mp}		  
   kHNLDcyNuMuMu   = 5,	 // N --> nu_{a}    \mu^{\mp} \mu^{\pm}. W and Z interfere
   kHNLDcyPiMu     = 6,  // N --> \pi^{\pm} \mu^{\mp}
   kHNLDcyPi0Pi0Nu = 7,	 // N --> \pi^{0}   \pi^{0} \nu (any kind)
   kHNLDcyPiPi0E   = 8,	 // N --> \pi^{\pm} \pi^{0}   e^{\mp}	  
   kHNLDcyPiPi0Mu  = 9,	 // N --> \pi^{\pm} \pi^{0} \mu^{\mp}	  
   kHNLDcyTEST     = 99, // N --> vv. Test only, not a valid FS.

 } HNLDecayMode_t;

  } // namespace HNL
} // namespace genie
#endif
