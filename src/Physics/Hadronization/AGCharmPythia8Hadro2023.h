//____________________________________________________________________________
/*!

\class    genie::AGCharmPythia8Hadro2023

\brief    Andreopoulos - Gallagher (AG) GENIE Charm Hadronization model.

          The model relies on empirical charm fragmentation and pT functions,
          as well as on experimentally-determined charm fractions, to produce
          the ID and 4-momentum of charmed hadron in charm production events.

          The remnant (non-charm) system is hadronised by a call to PYTHIA.

          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Hugh Gallagher <gallag@minos.phy.tufts.edu>
          Tufts University

\created  August 17, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _AGCHARM_PYTHIA8_HADRONIZATION_H_
#define _AGCHARM_PYTHIA8_HADRONIZATION_H_

#include "Physics/Hadronization/AGCharmPythiaBaseHadro2023.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Framework/Utils/Pythia8Singleton.h"
#endif

namespace genie {

class Spline;
class FragmentationFunctionI;

class AGCharmPythia8Hadro2023 : public AGCharmPythiaBaseHadro2023 {

public:
  AGCharmPythia8Hadro2023();
  AGCharmPythia8Hadro2023(string config);
  virtual ~AGCharmPythia8Hadro2023();


private:
  void           Initialize      (void)                                  const;
  bool           HadronizeRemnant(int qrkSyst1, int qrkSyst2, double WR, TLorentzVector p4R,
                                  unsigned int& rpos, TClonesArray * particle_list) const;

  void LoadConfig (void);

  void CopyOriginalDecayFlags     (void) const;
  void SetDesiredDecayFlags       (void) const;
  void RestoreOriginalDecayFlags  (void) const;

  // Original PYTHIA decay flags (stored so as to be restored after hadronization)
  mutable bool fOriDecayFlag_pi0; // pi^0
  mutable bool fOriDecayFlag_K0;  // K^0
  mutable bool fOriDecayFlag_K0b; // \bar{K^0}
  mutable bool fOriDecayFlag_L0;  // \Lambda^0
  mutable bool fOriDecayFlag_L0b; // \bar{\Lambda^0}
  mutable bool fOriDecayFlag_Dm;  // \Delta^-
  mutable bool fOriDecayFlag_D0;  // \Delta^0
  mutable bool fOriDecayFlag_Dp;  // \Delta^+
  mutable bool fOriDecayFlag_Dpp; // \Delta^++
  // Required PYTHIA decay flags set via Configure() [fixed]
  mutable bool fReqDecayFlag_pi0;         // pi^0
  mutable bool fReqDecayFlag_K0;          // K^0
  mutable bool fReqDecayFlag_K0b;         // \bar{K^0}
  mutable bool fReqDecayFlag_L0;          // \Lambda^0
  mutable bool fReqDecayFlag_L0b;         // \bar{\Lambda^0}
  mutable bool fReqDecayFlag_Dm;          // \Delta^-
  mutable bool fReqDecayFlag_D0;          // \Delta^0
  mutable bool fReqDecayFlag_Dp;          // \Delta^+
  mutable bool fReqDecayFlag_Dpp;         // \Delta^++

};

}         // genie namespace

#endif    // _AGCHARM_PYTHIA8_HADRONIZATION__H_
