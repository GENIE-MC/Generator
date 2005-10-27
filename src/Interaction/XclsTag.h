//____________________________________________________________________________
/*!

\class    genie::XclsTag

\brief    Contains minimal information for tagging exclusive processes.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  December 08, 2004

*/
//___________________________________________________________________________

#ifndef _FINAL_STATE_H_
#define _FINAL_STATE_H_

#include <iostream>
#include <string>

#include <TObject.h>

#include "BaryonResonance/BaryonResonance.h"

using std::ostream;
using std::string;

namespace genie {

class XclsTag : public TObject {

public:

  XclsTag();
  XclsTag(const XclsTag & xcls);
  virtual ~XclsTag();

  //-- Getting Exclusive Final State information

  bool IsCharmEvent       (void) const { return fIsCharmEvent;     }
  bool IsInclusiveCharm   (void) const;
  int  CharmHadronPDGCode (void) const { return fCharmedHadronPdg; }

  int  NProtons           (void) const { return fNProtons;  }
  int  NNeutrons          (void) const { return fNNeutrons; }
  int  NPi0               (void) const { return fNPi0;      }
  int  NPiPlus            (void) const { return fNPiPlus;   }
  int  NPiMinus           (void) const { return fNPiMinus;  }
  int  NNucleons          (void) const { return fNNeutrons + fNProtons;       }
  int  NPions             (void) const { return fNPi0 + fNPiPlus + fNPiMinus; }

  bool        KnownResonance (void) const { return (fResonance != kNoResonance); }
  Resonance_t Resonance      (void) const { return fResonance; }

  //-- Setting Exclusive Final State information

  void SetCharm       (int charm_pdgc = 0);
  void SetNPions      (int npi_plus, int npi_0, int npi_minus);
  void SetNNucleons   (int np, int nn);
  void SetNProtons    (int np) { fNProtons  = np; }
  void SetNNeutrons   (int nn) { fNNeutrons = nn; }
  void UnsetCharm     (void);
  void ResetNPions    (void);
  void ResetNNucleons (void);
  void SetResonance   (Resonance_t res);

  void   Copy     (const XclsTag & final_state);
  void   Print    (ostream & stream) const;
  string AsString (void) const;

  friend ostream & operator << (ostream& stream, const XclsTag & fin_state);

private:

  void Initialize(void);

  bool        fIsCharmEvent;     // true if we have charm production
  int         fCharmedHadronPdg; // charmed hadron pdg-code
  int         fNProtons;         // number of p's in the f/s hadronic system
  int         fNNeutrons;        // number of n's in the f/s hadronic system
  int         fNPi0;             // number of pi^0's in the f/s hadronic system
  int         fNPiPlus;          // number of pi^+'s in the f/s hadronic system
  int         fNPiMinus;         // number of pi^-'s in the f/s hadronic system
  Resonance_t fResonance;        // baryon resonance excited by probe

ClassDef(XclsTag,1)
};

}      // namespace

#endif // _FINAL_STATE_H_
