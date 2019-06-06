//____________________________________________________________________________
/*!

\class    genie::NucleonDecayMode

\brief    Enumeration of nucleon decay modes.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 10, 2011

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _NUCLEON_DECAY_MODE_H_
#define _NUCLEON_DECAY_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

 typedef enum ENucleonDecayMode {

   kNDNull = 0,
   // Antilepton + meson
   kNDN2eppi,           // m = 1: p --> e^{+} + \pi^{0}, n --> e^{+} + \pi^{-}
   kNDN2muppi,          // m = 2: p --> \mu^{+} + \pi^{0}, n --> \mu^{+} + \pi^{-}
   kNDN2nubarpi,        // m = 3: p --> \bar{\nu}} + \pi^{+}, n -->  \bar{\nu}} + \pi^{0} 
   kNDp2epeta,          // m = 4: p --> e^{+} + \eta
   kNDp2mupeta,         // m = 5: p --> \mu^{+} + \eta
   kNDn2nubareta,       // m = 6: n --> \bar{\nu}} + \eta
   kNDN2eprho,          // m = 7: p --> e^{+} + \rho^{0}, n --> e^{+} + \rho^{-}
   kNDN2muprho,         // m = 8: p --> \mu^{+} + \rho^{0}, n --> \mu^{+} + \rho^{-}
   kNDN2nubarrho,       // m = 9: p --> \bar{\nu}} + \rho^{+}, n --> \bar{\nu}} + \rho^{0} 
   kNDp2epomega,        // m = 10: p --> e^{+} + \omega
   kNDp2mupomega,       // m = 11: p --> \mu^{+} + \omega
   kNDn2nubaromega,     // m = 12: n --> \bar{\nu}} + \omega
   kNDN2epK,            // m = 13: p --> e^{+} + K^{0}, n --> e^{+} + K^{-}
   kNDp2epK0s,          // m = 14: p --> e^{+} + K^{0}_{short}
   kNDp2epK0l,          // m = 15: p --> e^{+} + K^{0}_{long}
   kNDN2mupK,           // m = 16: p --> \mu^{+} + K^{0}, n --> \mu^{+} + K^{-}
   kNDp2mupK0s,         // m = 17: p --> \mu^{+} + K^{0}_{short}
   kNDp2mupK0l,         // m = 18: p --> \mu^{+} + K^{0}_{long}
   kNDN2nubarK,         // m = 19: p --> \bar{\nu}} + K^{+}, n --> \bar{\nu}} + K^{0} 
   kNDn2nubarK0s,       // m = 20: n --> \bar{\nu}} + K^{0}_{short} 
   kNDp2epKstar0,       // m = 21: p --> e^{+} + K^{\star 0}
   kNDN2nubarKstar,     // m = 22: p --> \bar{\nu}} + K^{\star +}, n --> \bar{\nu}} + K^{\star 0}
   // Antilepton + mesons
   kNDp2eppippim,       // m = 23: p --> e^{+} + \pi^{+} + \pi^{-}
   kNDp2eppi0pi0,       // m = 24: p --> e^{+} + \pi^{0} + \pi^{0}
   kNDn2eppimpi0,       // m = 25: n --> e^{+} + \pi^{-} + \pi^{0}
   kNDp2muppippim,      // m = 26: p --> \mu^{+} + \pi^{+} + \pi^{-}
   kNDp2muppi0pi0,      // m = 27: p --> \mu^{+} + \pi^{0} + \pi^{0}
   kNDn2muppimpi0,      // m = 28: n --> \mu^{+} + \pi^{-} + \pi^{0}
   kNDn2epK0pim,        // m = 29: n --> e^{+} + K^{0} + \pi^{-}
   // Lepton + meson
   kNDn2empip,          // m = 30: n --> e^{-} + \pi^{+}
   kNDn2mumpip,         // m = 31: n --> \mu^{-} + \pi^{+}
   kNDn2emrhop,         // m = 32: n --> e^{-} + \rho^{+}
   kNDn2mumrhop,        // m = 33: n --> \mu^{-} + \rho^{+}
   kNDn2emKp,           // m = 34: n --> e^{-} + K^{+}
   kNDn2mumKp,          // m = 35: n --> \mu^{-} + K^{+}
   // Lepton + mesons
   kNDp2empippip,       // m = 36: p --> e^{-} + \pi^{+} + \pi^{+}
   kNDn2empippi0,       // m = 37: n --> e^{-} + \pi^{+} + \pi^{0}
   kNDp2mumpippip,      // m = 38: p --> \mu^{-} + \pi^{+} + \pi^{+}
   kNDn2mumpippi0,      // m = 39: n --> \mu^{-} + \pi^{+} + \pi^{0}
   kNDp2empipKp,        // m = 40: p --> e^{-} + \pi^{+} + K^{+}
   kNDp2mumpipKp,       // m = 41: p --> \mu^{-} + \pi^{+} + K^{+}
   // Antilepton + photon(s)
   kNDp2epgamma,        // m = 42: p --> e^{+} + \gamma
   kNDp2mupgamma,       // m = 43: p --> \mu^{+} + \gamma
   kNDn2nubargamma,     // m = 44: n -->  \bar{\nu}} + \gamma
   kNDp2epgammagamma,   // m = 45: p --> e^{+} + \gamma + \gamma
   kNDn2nubargammagamma,// m = 46: n -->  \bar{\nu}} + \gamma + \gamma
   // Three (or more) leptons
   kNDp2epepem = 49,    // m = 49: p --> e^{+} + e^{+} + e^{-}
   kNDp2epmupmum,       // m = 50: p --> e^{+} + \mu^{+} + \mu^{-}
   kNDp2epnubarnu,      // m = 51: p --> e^{+} + \bar{\nu}} + \nu
   kNDn2epemnubar,      // m = 52: n --> e^{+} + e^{-} + \bar{\nu}}
   kNDn2mupemnubar,     // m = 53: n --> \mu^{+} + e^{-} + \bar{\nu}}
   kNDn2mupmumnubar,    // m = 54: n --> \mu^{+} + \mu^{-} + \bar{\nu}}
   kNDp2mupepem,        // m = 55: p --> \mu^{+} + e^{+} + e^{-}
   kNDp2mupmupmum,      // m = 56: p --> \mu^{+} + \mu^{+} + \mu^{-}
   kNDp2mupnubarnu,     // m = 57: p --> \mu^{+} + \bar{\nu}} + \nu
   kNDp2emmupmup,       // m = 58: p --> e^{-} + \mu^{+} + \mu^{+}
   kNDn2threenus,       // m = 59: n --> \bar{\nu}} + \bar{\nu}} + \nu
   kNDn2fivenus         // m = 60: n --> \bar{\nu}} + \bar{\nu}} + \bar{\nu}} + \nu + \nu 

 } NucleonDecayMode_t;

}
#endif
