//____________________________________________________________________________
/*!

\namespace genie::intranuke

\brief     INTRANUKE utilities

\author    Jim Dobson <j.dobson07 \at imperial.ac.uk>
           Imperial College London

\created   Mar 03, 2009

\cpright   Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_UTILS_H_
#define _INTRANUKE_UTILS_H_

#include "HadronTransport/INukeHadroFates.h"

namespace genie {

class GHepRecord;

namespace utils {
namespace intranuke
{

  //! Reconstruct the INTRANUKE/hA model fate for the hadron at position i
  INukeFateHA_t ReconstructHadronFateHA  (GHepRecord * event, int i); 

}      // intranuke namespace
}      // utils     namespace
}      // genie     namespace

#endif // _INTRANUKE_UTILS_H_
