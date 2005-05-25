//____________________________________________________________________________
/*!

\class    genie::InitialState

\brief    Initial State information

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004

*/
//____________________________________________________________________________

#ifndef _INITIAL_STATE_H_
#define _INITIAL_STATE_H_

#include <iostream>
#include <string>

#include <TParticlePDG.h>
#include <TLorentzVector.h>

#include "Conventions/RefFrame.h"
#include "Interaction/Target.h"

using std::ostream;
using std::string;

namespace genie {

class InitialState {

public:

  InitialState();
  InitialState(int tgt_pdgc, int probe_pdgc);
  InitialState(int Z, int A, int probe_pdgc);
  InitialState(const Target & tgt, int probe_pdgc);
  InitialState(const InitialState & initial_state);
  virtual ~InitialState();

  const Target &   GetTarget       (void) const { return *fTarget;    }
  int              GetProbePDGCode (void) const { return  fProbePdgC; }
  Target *         GetTargetPtr    (void) const { return  fTarget;    }
  TParticlePDG *   GetProbe        (void) const;
  TLorentzVector * GetTargetP4 (RefFrame_t ref_frame = kRfLab) const;
  TLorentzVector * GetProbeP4  (RefFrame_t ref_frame = kRfStruckNucAtRest) const;
  double           GetProbeE   (RefFrame_t ref_frame) const;
   
  void SetProbePDGCode (int pdg_code);
  void SetTargetP4     (const TLorentzVector & P4); // in LAB-frame
  void SetProbeP4      (const TLorentzVector & P4); // in LAB-frame

  string AsString (void) const;
  void   Print    (ostream & stream) const;
  
  friend ostream & operator << (ostream& stream, const InitialState & init_state);
   
private:

  void Initialize (void);
  void Copy       (const InitialState & init_state);
  void Create     (int Z, int A, int probe_pdgc);

  int              fProbePdgC; // probe PDG code
  Target *         fTarget;    // nuclear target
  TLorentzVector * fProbeP4;   // probe 4-momentum in LAB-frame
  TLorentzVector * fTargetP4;  // nuclear target 4-momentum in LAB-frame
};

}      // namespace

#endif // _INITIAL_STATE_H_
