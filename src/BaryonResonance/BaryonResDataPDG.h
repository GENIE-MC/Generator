//____________________________________________________________________________
/*!

\class    genie::BaryonResDataPDG

\brief    Concrete implementation of the BaryonResDataSetI interface: Its
          configuration registry is loaded from an XML file with PDG baryon
          resonance data and they are served on request.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _BARYON_RES_DATA_PDG_H_
#define _BARYON_RES_DATA_PDG_H_

#include "BaryonResonance/BaryonResDataSetI.h"

namespace genie {

class BaryonResDataPDG : public BaryonResDataSetI {

public:

  BaryonResDataPDG();
  BaryonResDataPDG(string config);
  virtual ~BaryonResDataPDG();

  //-- BaryonResDataSetI interface implementation

  int    ResonanceIndex    (Resonance_t res) const;
  int    OrbitalAngularMom (Resonance_t res) const;
  bool   IsDeltaResonance  (Resonance_t res) const;
  bool   IsNResonance      (Resonance_t res) const;
  double Mass              (Resonance_t res) const;
  double Width             (Resonance_t res) const;
  double BreitWignerNorm   (Resonance_t res) const;

};

}         // genie namespace 

#endif    // _BARYON_RES_DATA_PDG_H_
