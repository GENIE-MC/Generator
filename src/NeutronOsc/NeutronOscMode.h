//____________________________________________________________________________
/*!

\class    genie::NeutronOscMode

\brief    Enumeration of neutron oscillation annihilation modes.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 10, 2011

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _NEUTRON_OSC_MODE_H_
#define _NEUTRON_OSC_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

 typedef enum ENeutronOscMode {

     // i just replaced all the nucleon decay modes with nnbar modes -j

     kNONull     = -1,
     kNORandom,            // Will select a random decay mode -j
     kNOpto1pip1pi0,       // p + nbar --> \pi^{+} + \pi^{0}
     kNOpto1pip2pi0,       // p + nbar --> \pi^{+} + 2\pi^{0}
     kNOpto1pip3pi0,       // p + nbar --> \pi^{+} + 3\pi^{0}
     kNOpto2pip1pim1pi0,   // p + nbar --> 2\pi^{+} + \pi^{-} + \pi^{0}
     kNOpto2pip1pim2pi0,   // p + nbar --> 2\pi^{+} + \pi^{-} + 2\pi^{0}
     kNOpto2pip1pim2o,     // p + nbar --> 2\pi^{+} + \pi^{-} + 2\omega^{0}
     kNOpto3pip2pim1pi0,   // p + nbar --> 3\pi^{+} + 2\pi^{-} + \pi^{0}
     kNOnto1pip1pim,       // n + nbar --> \pi^{+} + \pi^{-}
     kNOnto2pi0,           // n + nbar --> 2\pi^{0}
     kNOnto1pip1pim1pi0,   // n + nbar --> \pi^{+} + \pi^{-} + \pi^{0}
     kNOnto1pip1pim2pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 2\pi^{0}
     kNOnto1pip1pim3pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 3\pi^{0}
     kNOnto2pip2pim,       // n + nbar --> 2\pi^{+} + 2\pi^{-}
     kNOnto2pip2pim1pi0,   // n + nbar --> 2\pi^{+} + 2\pi^{-} + \pi^{0}
     kNOnto1pip1pim1o,     // n + nbar --> \pi^{+} + \pi^{-} + \omega^{0}
     kNOnto2pip2pim2pi0    // n + nbar --> 2\pi^{+} + 2\pi^{-} + 2\pi^{0}

 } NeutronOscMode_t;

}
#endif
