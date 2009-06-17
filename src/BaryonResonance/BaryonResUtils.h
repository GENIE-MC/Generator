//____________________________________________________________________________
/*!

\namespace genie::utils::res

\brief     Baryon Resonance utilities.

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   November 25, 2004

\cpright   Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_UTILS_H_
#define _BARYON_RESONANCE_UTILS_H_

#include <string>

#include "BaryonResonance/BaryonResonance.h"

using std::string;

namespace genie {

class Interaction;

namespace utils {
namespace res {

  const char* AsString          (Resonance_t res);
  double      Mass              (Resonance_t res);
  Resonance_t FromString        (const char * res);
  Resonance_t FromPdgCode       (int pdgc);
  int         PdgCode           (Resonance_t res, int Q);
  bool        IsBaryonResonance (int pdgc);
  bool        IsDelta           (Resonance_t res);
  bool        IsN               (Resonance_t res);
  int         ResonanceCharge   (const Interaction * interaction);

}        // res   namespace
}        // utils namespace
}        // genie namespace

#endif   // _BARYON_RESONANCE_UTILS_H_
