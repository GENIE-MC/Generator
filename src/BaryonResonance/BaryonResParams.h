//____________________________________________________________________________
/*!

\class    genie::BaryonResParams

\brief    Encapsulates baryon resonance parameters. Uses the Strategy pattern
          to retrieve data from a concrete BaryonResDataSetI.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _BARYON_RES_PARAMS_H_
#define _BARYON_RES_PARAMS_H_

#include <iostream>

#include "BaryonResonance/BaryonResonance.h"

using std::ostream;

namespace genie {

class BaryonResDataSetI;

class BaryonResParams
{
public:

  BaryonResParams();
  BaryonResParams(const BaryonResParams & res_params);
  virtual ~BaryonResParams();

  void   SetDataSet   (const BaryonResDataSetI * data_set);
  void   RetrieveData (Resonance_t resonance);
    
  int    ResonanceIndex    (void) const { return fResonanceIndex;    }
  int    OrbitalAngularMom (void) const { return fOrbitalAngularMom; }
  bool   IsDeltaResonance  (void) const { return fIsDeltaResonance;  }
  bool   IsNResonance      (void) const { return fIsNResonance;      }
  double Mass              (void) const { return fMass;              }
  double Width             (void) const { return fWidth;             }
  double BreitWignerNorm   (void) const { return fBreitWignerNorm;   }

  void   Print(ostream & stream) const;
  
  friend ostream & operator<<(ostream & stream, const BaryonResParams & res_params);

private:

  void   InitParams();
  
  int    fResonanceIndex;     ///< resonance index (quantum number n)
  int    fOrbitalAngularMom;  ///< resonannce orbital angular momentum L
  bool   fIsDeltaResonance;   ///< true if it is a Delta resonance
  bool   fIsNResonance;       ///< true if it is a N resonance
  double fMass;               ///< resonance Breit-Wigner mass
  double fWidth;              ///< resonance Breit-Wigner width
  double fBreitWignerNorm;    ///< resonance Breit-Wigner normalization

  const BaryonResDataSetI * fDataSet;  
};

}       // genie namepace

#endif  // _BARYON_RES_PARAMS_H_
