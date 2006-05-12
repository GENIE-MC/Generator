//____________________________________________________________________________
/*!

\class    genie::InitialState

\brief    Initial State information

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INITIAL_STATE_H_
#define _INITIAL_STATE_H_

#include <iostream>
#include <string>

#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TObject.h>

#include "Conventions/RefFrame.h"
#include "Interaction/Target.h"

using std::ostream;
using std::string;

namespace genie {

class InitialState : public TObject {

public:

  InitialState();
  InitialState(int tgt_pdgc, int probe_pdgc);
  InitialState(int Z, int A, int probe_pdgc);
  InitialState(const Target & tgt, int probe_pdgc);
  InitialState(const InitialState & initial_state);
  ~InitialState();

  TParticlePDG *   GetProbe         (void) const;
  const Target &   GetTarget        (void) const { return *fTarget;    }
  Target *         GetTargetPtr     (void) const { return  fTarget;    }
  int              GetProbePDGCode  (void) const { return  fProbePdgC; }
  int              GetTargetPDGCode (void) const;

  TLorentzVector * GetTargetP4 (RefFrame_t ref_frame = kRfLab) const;
  TLorentzVector * GetProbeP4  (RefFrame_t ref_frame = kRfStruckNucAtRest) const;
  double           GetProbeE   (RefFrame_t ref_frame) const;

  void SetPDGCodes      (int tgt_pdgc, int probe_pdgc);
  void SetProbePDGCode  (int pdg_code);
  void SetTargetPDGCode (int pdg_code);
  void SetTargetP4      (const TLorentzVector & P4); // in LAB-frame
  void SetProbeP4       (const TLorentzVector & P4); // in LAB-frame

  //! Copy, reset, compare, print itself and build string code
  void   Reset    (void);
  void   Copy     (const InitialState & init_state);
  bool   Compare  (const InitialState & init_state) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const InitialState & i) const;
  InitialState &   operator =  (const InitialState & i);
  friend ostream & operator << (ostream & stream, const InitialState & i);

private:

  //! Methods for InitialState initialization and clean up
  void Init       (void);
  void Init       (int target_pdgc, int probe_pdgc);
  void CleanUp    (void);

  //! Private data members
  int              fProbePdgC; ///< probe PDG code
  Target *         fTarget;    ///< nuclear target
  TLorentzVector * fProbeP4;   ///< probe 4-momentum in LAB-frame
  TLorentzVector * fTargetP4;  ///< nuclear target 4-momentum in LAB-frame

ClassDef(InitialState,1)
};

}      // namespace

#endif // _INITIAL_STATE_H_
