//____________________________________________________________________________
/*!

\class    genie::KPhaseSpaceCuts

\brief    Singleton service for configurable kinematic phase-space cuts.

\author   GENIE Collaboration

\created  May 21, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _KINEMATIC_PHASE_SPACE_CUTS_H_
#define _KINEMATIC_PHASE_SPACE_CUTS_H_

#include <string>

using std::string;

namespace genie {

class Interaction;

class KPhaseSpaceCuts {

public:
  static KPhaseSpaceCuts * Instance(void);

  void   SetQ2MinOverride(double q2min);
  bool   HasQ2MinCut(const Interaction * interaction) const;
  double Q2MinCut(const Interaction * interaction, double default_q2min) const;
  bool   HasSplineQ2MinCut(void) const;
  double SplineQ2MinCut(void) const;
  string SplineQ2MinCutSource(void) const;
  bool   HasSplineQ2MinCutForProbe(int probe_pdg) const;
  double SplineQ2MinCutForProbe(int probe_pdg) const;
  string SplineQ2MinCutSourceForProbe(int probe_pdg) const;

private:
  KPhaseSpaceCuts();
  KPhaseSpaceCuts(const KPhaseSpaceCuts & cuts);
  virtual ~KPhaseSpaceCuts();

  void   LoadConfig(void) const;
  double ValidateQ2Min(const char * name, double q2min) const;
  double RaiseDefault(double default_q2min, double cut_q2min) const;
  void   FailMissingEMQ2Min(void) const;

  static KPhaseSpaceCuts * fInstance;

  mutable bool   fLoaded;
  mutable bool   fHasEMQ2Min;
  mutable bool   fHasWeakQ2Min;
  mutable double fEMQ2Min;
  mutable double fWeakQ2Min;

  bool   fHasQ2MinOverride;
  double fQ2MinOverride;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (KPhaseSpaceCuts::fInstance != 0) {
            delete KPhaseSpaceCuts::fInstance;
            KPhaseSpaceCuts::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _KINEMATIC_PHASE_SPACE_CUTS_H_
