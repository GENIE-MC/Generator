//____________________________________________________________________________
/*!

  \class    genie::NNBarOscMode

  \brief    Enumeration of neutron oscillation annihilation modes.

  \author   Jeremy Hewes, Georgia Karagiorgi
  University of Manchester

  \created  November, 2016

  \author   Linyan Wan
  Fermilab

  \created  April, 2024


  \cpright  Copyright (c) 2003-2023, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org         
  */
//____________________________________________________________________________

#ifndef _N_NBAR_OSC_MODE_H_
#define _N_NBAR_OSC_MODE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

  typedef enum ENNBarOscMode {

	 // i just replaced all the nucleon decay modes with nnbar modes -j

	 kNONull     = -1,
	 kNORandom,            // Will select a random decay mode -j
	 kNOpto1pip1pi0,       // p + nbar --> \pi^{+} + \pi^{0}
			 kNOpto1pip2pi0,       // p + nbar --> \pi^{+} + 2\pi^{0}
			 kNOpto1pip3pi0,       // p + nbar --> \pi^{+} + 3\pi^{0}
			 kNOpto1pip4pi0,       // p + nbar --> \pi^{+} + 4\pi^{0}
			 kNOpto2pip1pim,   		// p + nbar --> 2\pi^{+} + \pi^{-}
			 kNOpto2pip1pim1pi0,   // p + nbar --> 2\pi^{+} + \pi^{-} + \pi^{0}
			 kNOpto2pip1pim2pi0,   // p + nbar --> 2\pi^{+} + \pi^{-} + 2\pi^{0}
			 kNOpto2pip1pim3pi0,   // p + nbar --> 2\pi^{+} + \pi^{-} + 3\pi^{0}
			 kNOpto3pip2pim,		   // p + nbar --> 3\pi^{+} + 2\pi^{-}
			 kNOpto3pip2pim1pi0,   // p + nbar --> 3\pi^{+} + 2\pi^{-} + \pi^{0}
			 kNOpto1pip1pi01o,     // p + nbar --> 1\pi^{+} + \pi^{0} + 1\omega^{0}
			 kNOpto2pip1pim1o,     // p + nbar --> 2\pi^{+} + \pi^{-} + 1\omega^{0}

			 kNOnto2pi0,           // n + nbar --> 2\pi^{0}
			 kNOnto3pi0,           // n + nbar --> 3\pi^{0}
			 kNOnto4pi0,           // n + nbar --> 4\pi^{0}
			 kNOnto5pi0,           // n + nbar --> 5\pi^{0}
			 kNOnto7pi0,           // n + nbar --> 7\pi^{0}
			 kNOnto1pip1pim,       // n + nbar --> \pi^{+} + \pi^{-}
			 kNOnto1pip1pim1pi0,   // n + nbar --> \pi^{+} + \pi^{-} + \pi^{0}
			 kNOnto1pip1pim2pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 2\pi^{0}
			 kNOnto1pip1pim3pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 3\pi^{0}
			 kNOnto1pip1pim4pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 4\pi^{0}
			 kNOnto1pip1pim5pi0,   // n + nbar --> \pi^{+} + \pi^{-} + 5\pi^{0}
			 kNOnto2pip2pim,       // n + nbar --> 2\pi^{+} + 2\pi^{-}
			 kNOnto2pip2pim1pi0,   // n + nbar --> 2\pi^{+} + 2\pi^{-} + \pi^{0}
			 kNOnto2pip2pim2pi0,   // n + nbar --> 2\pi^{+} + 2\pi^{-} + 2\pi^{0}
			 kNOnto2pip2pim3pi0,   // n + nbar --> 2\pi^{+} + 2\pi^{-} + 3\pi^{0}
			 kNOnto3pip3pim,   		// n + nbar --> 3\pi^{+} + 3\pi^{-}
			 kNOnto3pip3pim1pi0,   // n + nbar --> 3\pi^{+} + 3\pi^{-} + 1\pi^{0}
			 kNOnto1rho1pi0,			// n + nbar --> \rho^{0} + \pi^{0}
			 kNOnto1rhopm1pimp,		// n + nbar --> \rho^{+-} + \pi^{-+}
			 kNOnto2o,     			// n + nbar --> 2\omega^{0}
			 kNOnto1rho1o,     		// n + nbar --> \rho^{0} + \omega^{0}
			 kNOnto2pi01o,     		// n + nbar --> 2\pi^{0} + \omega^{0}
			 kNOnto1pip1pim1o,     // n + nbar --> \pi^{+} + \pi^{-} + \omega^{0}
			 kNOnto1eta1o,     		// n + nbar --> \eta^{0} + \omega^{0}
			 kNOnto1pip1pim1eta,   // n + nbar --> \pi^{+} + \pi^{-} + \eta^{0}
			 kNOnto1Kp1Km,     		// n + nbar --> \K^{+} + \K^{-}
			 kNOnto1Kp1Km1o,     	// n + nbar --> K^{+} + K^{-} + \omega^{0}
			 kNOnto2Ks1o,     		// n + nbar --> 2K^{s} + \omega^{0}
			 kNOnto1pipm1Kmp1Kl,   // n + nbar --> \pi^{+-} + K^{-+} + K^{l}
			 kNOnto1pipm1Kmp1Ks,   // n + nbar --> \pi^{+-} + K^{-+} + K^{s}
			 kNOnto1pipm1Kmp1pi0K0,// n + nbar --> \pi^{+-} + K^{-+} + \pi^{0} + K^{0}
			 kNOnto1pip1pim1Ks1Kl   // n + nbar --> \pi^{+} + \pi^{-} + K^{s} + K^{l}

  } NNBarOscMode_t;

}
#endif
