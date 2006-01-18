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
  ~XclsTag();

  //! Getting Exclusive Intermediate and/or Final State information
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
  bool KnownResonance     (void) const { return (fResonance != kNoResonance); }
  Resonance_t Resonance   (void) const { return fResonance; }

  //! Setting Exclusive Final State information
  void SetCharm       (int charm_pdgc = 0);
  void SetNPions      (int npi_plus, int npi_0, int npi_minus);
  void SetNNucleons   (int np, int nn);
  void SetNProtons    (int np) { fNProtons  = np; }
  void SetNNeutrons   (int nn) { fNNeutrons = nn; }
  void UnsetCharm     (void);
  void ResetNPions    (void);
  void ResetNNucleons (void);
  void SetResonance   (Resonance_t res);

  //! Copy, reset, print itself and build string code
  void   Reset    (void);
  void   Copy     (const XclsTag & final_state);
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  XclsTag &        operator =  (const XclsTag & xcls);
  friend ostream & operator << (ostream& stream, const XclsTag & xcls);

private:

  //! Private data members
  bool        fIsCharmEvent;     ///< true if we have charm production
  int         fCharmedHadronPdg; ///< charmed hadron pdg-code
  int         fNProtons;         ///< # of p's in the f/s hadronic system
  int         fNNeutrons;        ///< # of n's in the f/s hadronic system
  int         fNPi0;             ///< # of pi^0's in the f/s hadronic system
  int         fNPiPlus;          ///< # of pi^+'s in the f/s hadronic system
  int         fNPiMinus;         ///< # of pi^-'s in the f/s hadronic system
  Resonance_t fResonance;        ///< baryon resonance excited by probe

ClassDef(XclsTag,1)
};

}      // namespace

#endif // _FINAL_STATE_H_
