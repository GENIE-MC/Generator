//____________________________________________________________________________
/*!

\class    genie::BaryonResDataPDG

\brief    Concrete implementation of the BaryonResDataSetI interface: Its
          configuration registry is loaded from an XML file with PDG baryon
          resonance data and they are served on request.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BARYON_RES_DATA_PDG_H_
#define _BARYON_RES_DATA_PDG_H_

#include <map>

#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResList.h"

using std::map;

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

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadResonanceData(void);

  BaryonResList fResList;

  map<Resonance_t, int>    fResIdx;
  map<Resonance_t, int>    fResL;
  map<Resonance_t, bool>   fIsD;
  map<Resonance_t, bool>   fIsN;
  map<Resonance_t, double> fResMass;
  map<Resonance_t, double> fResWidth;
  map<Resonance_t, double> fResNorm;
};

}         // genie namespace

#endif    // _BARYON_RES_DATA_PDG_H_
