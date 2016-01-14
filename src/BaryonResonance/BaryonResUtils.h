//____________________________________________________________________________
/*!

\namespace genie::utils::res

\brief     Baryon Resonance utilities.

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

\created   November 25, 2004

\cpright   Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

namespace utils {
namespace res {

  const char* AsString   (Resonance_t res);  ///< resonance id -> string
  Resonance_t FromString (const char * res); ///< string -> resonance id

  int         PdgCode     (Resonance_t res, int Q); ///< (resonance id, charge) -> PDG code
  Resonance_t FromPdgCode (int pdgc);               ///< PDG code -> resonance id

  bool        IsBaryonResonance (int pdgc);         ///< is input a baryon resonance?
  bool        IsDelta           (Resonance_t res);  ///< is it a Delta resonance?
  bool        IsN               (Resonance_t res);  ///< is it an N resonance?
  double      Mass              (Resonance_t res);  ///< resonance mass (GeV)
  double      Width             (Resonance_t res);  ///< resonance width (GeV)
  double      BWNorm            (Resonance_t res);  ///< breit-wigner normalization factor
  int         OrbitalAngularMom (Resonance_t res);  ///< orbital angular momentum
  int         ResonanceIndex    (Resonance_t res);  ///< resonance idx, quark model / SU(6)

}        // res   namespace
}        // utils namespace
}        // genie namespace

#endif   // _BARYON_RESONANCE_UTILS_H_
