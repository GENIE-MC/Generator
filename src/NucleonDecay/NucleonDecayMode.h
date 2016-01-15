//____________________________________________________________________________
/*!

\class    genie::NucleonDecayMode

\brief    Enumeration of nucleon decay modes.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 10, 2011

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

     kNDNull     = -1,
     kNDp2epi0,        // p --> e^{+}   + \pi^{0}
     kNDp2mupi0,       // p --> \mu^{+} + \pi^{0}
     kNDp2eeta0,       // p --> e^{+}   + \eta^{0}
     kNDp2mueta0,      // p --> \mu^{+} + \eta^{0}
     kNDp2erho0,       // p --> e^{+}   + \rho^{0}
     kNDp2murho0,      // p --> \mu^{+} + \rho^{0}
     kNDp2eomega0,     // p --> e^{+}   + \omega^{0}
     kNDp2muomega0,    // p --> \mu^{+} + \omega^{0}
     kNDn2epim,        // n --> e^{+}   + \pi^{-}
     kNDn2mupim,       // n --> \mu^{+} + \pi^{-}
     kNDp2nubarKp      // p --> \bar{\nu}} + K^{+}

 } NucleonDecayMode_t;

}
#endif
