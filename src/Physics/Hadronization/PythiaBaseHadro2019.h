//____________________________________________________________________________
/*!

\class    genie::PythiaBaseHadro2019

\brief    Base class for the Pythia (6 and 8) hadronization modules in GENIE.
          In particular, the base class provides common checks and basic
          assignments of quark/diquark codes for a no frills interface to
          Pythia hadronization routines.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PYTHIA_HADRONIZATION_BASE_H_
#define _PYTHIA_HADRONIZATION_BASE_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Interaction/Interaction.h"

namespace genie {

class PythiaBaseHadro2019 : public EventRecordVisitorI {

protected:
  PythiaBaseHadro2019();
  PythiaBaseHadro2019(string name);
  PythiaBaseHadro2019(string name, string config);
  virtual ~PythiaBaseHadro2019();

  virtual void ProcessEventRecord         (GHepRecord* event)     const;
  virtual void MakeQuarkDiquarkAssignments(const Interaction* in) const;
  virtual bool AssertValidity             (const Interaction* in) const;
  virtual void Initialize                 (void);
  virtual void LoadConfig                 (void);

  virtual void CopyOriginalDecayFlags     (void) const = 0;
  virtual void SetDesiredDecayFlags       (void) const = 0;
  virtual void RestoreOriginalDecayFlags  (void) const = 0;

  virtual bool Hadronize (GHepRecord* event) const = 0;

  // PDG codes assigned on an event-by-event basis and driving
  // PYTHIA6/8 hadronization routines
  mutable int fLeadingQuark;
  mutable int fRemnantDiquark;

  // PYTHIA physics configuration parameters used
  double fSSBarSuppression;       ///< ssbar suppression
  double fGaussianPt2;            ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail;     ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;       ///< remaining E cutoff stopping fragmentation
  double fDiQuarkSuppression;     ///< di-quark suppression parameter
  double fLightVMesonSuppression; ///< light vector meson suppression
  double fSVMesonSuppression;     ///< strange vector meson suppression
  double fLunda;                  ///< Lund a parameter
  double fLundb;                  ///< Lund b parameter
  double fLundaDiq;               ///< adjustment of Lund a for di-quark

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
  bool fReqDecayFlag_pi0;         // pi^0
  bool fReqDecayFlag_K0;          // K^0
  bool fReqDecayFlag_K0b;         // \bar{K^0}
  bool fReqDecayFlag_L0;          // \Lambda^0
  bool fReqDecayFlag_L0b;         // \bar{\Lambda^0}
  bool fReqDecayFlag_Dm;          // \Delta^-
  bool fReqDecayFlag_D0;          // \Delta^0
  bool fReqDecayFlag_Dp;          // \Delta^+
  bool fReqDecayFlag_Dpp;         // \Delta^++
};

}         // genie namespace

#endif    // _PYTHIA_HADRONIZATION_BASE_H_
