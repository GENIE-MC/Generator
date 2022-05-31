//____________________________________________________________________________
/*!

\class    genie::NHL::NHLDecayMode

\brief    Enumeration of NHL decay modes.

\author

\created  November 10, 2011

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NHL_DECAY_MODE_H_
#define _NHL_DECAY_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
  namespace NHL {

 typedef enum ENHLDecayMode {

   kNHLDcyNull     = -1, // dummy
   kNHLDcyNuNuNu   = 0,	 // N --> 3 nus. Summed over all flavours 
   kNHLDcyNuEE     = 1,	 // N --> nu_{a}      e^{\mp}   e^{\pm}. W and Z interfere
   kNHLDcyNuMuE    = 2,	 // N --> nu_{\mu/e}  e^{\mp} \mu^{\pm}. Only W. Summed over nue and numu
   kNHLDcyPi0Nu    = 3,	 // N --> \pi^{0}   \nu (any kind)	  
   kNHLDcyPiE      = 4,  // N --> \pi^{\pm}   e^{\mp}		  
   kNHLDcyNuMuMu   = 5,	 // N --> nu_{a}    \mu^{\mp} \mu^{\pm}. W and Z interfere
   kNHLDcyPiMu     = 6,  // N --> \pi^{\pm} \mu^{\mp}
   kNHLDcyPi0Pi0Nu = 7,	 // N --> \pi^{0}   \pi^{0} \nu (any kind)
   kNHLDcyPiPi0E   = 8,	 // N --> \pi^{\pm} \pi^{0}   e^{\mp}	  
   kNHLDcyPiPi0Mu  = 9,	 // N --> \pi^{\pm} \pi^{0} \mu^{\mp}	  
   kNHLDcyTEST     = 99, // N --> vv. Test only, not a valid FS.

 } NHLDecayMode_t;

  } // namespace NHL
} // namespace genie
#endif
