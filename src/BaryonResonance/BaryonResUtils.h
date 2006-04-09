//____________________________________________________________________________
/*!

\namespace genie::utils::res

\brief     Baryon Resonance utilities.

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   November 25, 2004

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

  char *      AsString          (Resonance_t res);
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
