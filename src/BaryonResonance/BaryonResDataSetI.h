//____________________________________________________________________________
/*!

\class    genie::BaryonResDataSetI

\brief    Pure abstract base class. Defines the BaryonResDataSetI interface
          to be implemented by any algorithmic class computing (or merely
          retrieving) baryon resonance data.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
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
