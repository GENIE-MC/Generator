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

#include <TParticlePDG.h>
#include <TLorentzVector.h>

#include "Conventions/RefFrame.h"
#include "Interaction/Target.h"

using std::ostream;

namespace genie {

class InitialState {

public:

  InitialState();
  InitialState(const Target & tgt, int probe_pdgc);
  InitialState(const InitialState & initial_state);
  virtual ~InitialState();

  const Target & GetTarget       (void) const { return *fTarget;    }
  int            GetProbePDGCode (void) const { return  fProbePdgC; }
  Target *       GetTargetPtr    (void) const { return  fTarget;    }
  TParticlePDG * GetProbe        (void) const;
  
  TLorentzVector * GetTargetP4 (RefFrame_t ref_frame = kRfLab) const;
  TLorentzVector * GetProbeP4  (RefFrame_t ref_frame = kRfStruckNucAtRest) const;
  
  void SetProbePDGCode (int pdg_code) { fProbePdgC = pdg_code; }

  void SetTargetP4     (const TLorentzVector & P4); // in LAB-frame
  void SetProbeP4      (const TLorentzVector & P4); // in LAB-frame

  void Copy  (const InitialState & init_state);
  void Print (ostream & stream) const;
  
  friend ostream & operator << (ostream& stream, const InitialState & init_state);
   
private:

  void Initialize(void);

  int              fProbePdgC; // probe PDG code
  Target *         fTarget;    // nuclear target
  TLorentzVector * fProbeP4;   // probe 4-momentum in LAB-frame
  TLorentzVector * fTargetP4;  // nuclear target 4-momentum in LAB-frame

};

}      // namespace

#endif // _INITIAL_STATE_H_
