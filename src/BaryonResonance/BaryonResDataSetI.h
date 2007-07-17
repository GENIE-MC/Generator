//____________________________________________________________________________
/*!

\class    genie::BaryonResDataSetI

\brief    Pure abstract base class. Defines the BaryonResDataSetI interface
          to be implemented by any algorithmic class computing (or merely
          retrieving) baryon resonance data.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _BARYON_RES_DATA_SET_I_H_
#define _BARYON_RES_DATA_SET_I_H_

#include "Algorithm/Algorithm.h"
#include "BaryonResonance/BaryonResonance.h"

namespace genie {

class BaryonResDataSetI : public Algorithm {

public:
  virtual ~BaryonResDataSetI();

  //-- define the BaryonResDataSetI interface
  
  virtual int    ResonanceIndex    (Resonance_t res) const = 0;
  virtual int    OrbitalAngularMom (Resonance_t res) const = 0;
  virtual bool   IsDeltaResonance  (Resonance_t res) const = 0;
  virtual bool   IsNResonance      (Resonance_t res) const = 0;
  virtual double Mass              (Resonance_t res) const = 0;
  virtual double Width             (Resonance_t res) const = 0;
  virtual double BreitWignerNorm   (Resonance_t res) const = 0;

protected:
  BaryonResDataSetI();
  BaryonResDataSetI(string name);
  BaryonResDataSetI(string name, string config);
};

}         // genie namespace 
#endif    // _BARYON_RES_DATA_SET_I_H_
