//____________________________________________________________________________
/*!

\class    genie::InitialState

\brief    Initial State information

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

          Other minor changes / additions and fixes were installed by: 
          Andy Furmanski (Univ. of Manchester)
          Joe Johnston (Univ of Pittsburgh)

\created  May 02, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
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

#include "Framework/Conventions/RefFrame.h"
#include "Framework/Interaction/Target.h"

using std::ostream;
using std::string;

class TRootIOCtor;

namespace genie {

class InitialState;
ostream & operator << (ostream & stream, const InitialState & i); 

class InitialState : public TObject {

public:
  using TObject::Print;   // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TObject::Copy;    // 
  using TObject::Compare; //

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
  double           CMEnergy   () const; ///< centre-of-mass energy (sqrt s)

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
  bool IsDMP    (void) const; ///< is dark matter   + proton?
  bool IsDMN    (void) const; ///< is dark matter   + neutron?

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
