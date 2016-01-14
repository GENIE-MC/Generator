//____________________________________________________________________________
/*!

\class    genie::InitialState

\brief    Initial State information

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 02, 2004

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

class TRootIOCtor;

namespace genie {

class InitialState : public TObject {

public:
  InitialState();
  InitialState(int tgt_pdgc, int probe_pdgc);
  InitialState(int Z, int A, int probe_pdgc);
  InitialState(const Target & tgt, int probe_pdgc);
  InitialState(const InitialState & initial_state);
  InitialState(TRootIOCtor*);
 ~InitialState();

  TParticlePDG *   Probe      (void) const;
  int              ProbePdg   (void) const { return fProbePdg; }
  int              TgtPdg     (void) const;
  const Target &   Tgt        (void) const { return *fTgt; }
  Target *         TgtPtr     (void) const { return  fTgt; }
  TLorentzVector * GetTgtP4   (RefFrame_t rf = kRfLab) const;
  TLorentzVector * GetProbeP4 (RefFrame_t rf = kRfHitNucRest) const;
  double           ProbeE     (RefFrame_t rf) const;

  void SetPdgs     (int tgt_pdgc, int probe_pdgc);
  void SetProbePdg (int pdg_code);
  void SetTgtPdg   (int pdg_code);
  void SetTgtP4    (const TLorentzVector & P4); // in LAB-frame
  void SetProbeP4  (const TLorentzVector & P4); // in LAB-frame
  void SetProbeE   (double E);                  // in LAB-frame (0,0,E,E)

  bool IsNuP    (void) const; ///< is neutrino      + proton?
  bool IsNuN    (void) const; ///< is neutrino      + neutron?
  bool IsNuBarP (void) const; ///< is anti-neutrino + proton?
  bool IsNuBarN (void) const; ///< is anti-neutrino + neutron?

  //-- Copy, reset, compare, print itself and build string code
  void   Reset    (void);
  void   Copy     (const InitialState & init_state);
  bool   Compare  (const InitialState & init_state) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  //-- Overloaded operators
  bool             operator == (const InitialState & i) const;             ///< equal?
  InitialState &   operator =  (const InitialState & i);                   ///< copy
  friend ostream & operator << (ostream & stream, const InitialState & i); ///< print

private:

  //-- Methods for InitialState initialization and clean up
  void Init       (void);
  void Init       (int target_pdgc, int probe_pdgc);
  void CleanUp    (void);

  //-- Private data members
  int              fProbePdg; ///< probe PDG code
  Target *         fTgt;      ///< nuclear target
  TLorentzVector * fProbeP4;  ///< probe 4-momentum in LAB-frame
  TLorentzVector * fTgtP4;    ///< nuclear target 4-momentum in LAB-frame

ClassDef(InitialState,1)
};

}      // namespace

#endif // _INITIAL_STATE_H_
