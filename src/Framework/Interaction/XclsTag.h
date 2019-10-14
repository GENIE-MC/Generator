//____________________________________________________________________________
/*!

\class    genie::XclsTag

\brief    Contains minimal information for tagging exclusive processes.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
          Marco Roda <mroda@liverpool.ac.uk>
          University of Liverpool

\created  December 08, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//___________________________________________________________________________

#ifndef _FINAL_STATE_H_
#define _FINAL_STATE_H_

#include <iostream>
#include <string>

#include <TObject.h>

#include "Framework/ParticleData/BaryonResonance.h"

using std::ostream;
using std::string;

namespace genie {

class XclsTag;
ostream & operator << (ostream& stream, const XclsTag & xcls); 

class XclsTag : public TObject {
  
public:
  using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TObject::Copy;

  XclsTag();
  XclsTag(const XclsTag & xcls);
 ~XclsTag();

  // Getting exclusive intermediate and/or final state information
  bool IsCharmEvent       (void) const { return fIsCharmEvent;     }
  bool IsInclusiveCharm   (void) const;
  int  CharmHadronPdg     (void) const { return fCharmedHadronPdg; }
  bool IsStrangeEvent     (void) const { return fIsStrangeEvent;   }
  bool IsInclusiveStrange (void) const;
  int  StrangeHadronPdg   (void) const {return fStrangeHadronPdg;  }
  int  NProtons           (void) const { return fNProtons;  }
  int  NNeutrons          (void) const { return fNNeutrons; }
  int  NPi0               (void) const { return fNPi0;      }
  int  NPiPlus            (void) const { return fNPiPlus;   }
  int  NPiMinus           (void) const { return fNPiMinus;  }
  int  NNucleons          (void) const { return fNNeutrons + fNProtons;       }
  int  NPions             (void) const { return fNPi0 + fNPiPlus + fNPiMinus; }
  int  NRhos              (void) const { return fNRho0 + fNRhoPlus + fNRhoMinus; }
  int  NSingleGammas      (void) const { return fNSingleGammas; }
  int  NRho0              (void) const { return fNRho0;      }
  int  NRhoPlus           (void) const { return fNRhoPlus;   }
  int  NRhoMinus          (void) const { return fNRhoMinus;  }
  bool KnownResonance     (void) const { return (fResonance != kNoResonance); }
  Resonance_t Resonance   (void) const { return fResonance; }
  int  DecayMode          (void) const { return fDecayMode; }

  // Ssetting exclusive final state information
  void SetCharm           (int charm_pdgc = 0);
  void SetStrange         (int strange_pdgc = 0);
  void SetNPions          (int npi_plus, int npi_0, int npi_minus);
  void SetNNucleons       (int np, int nn);
  void SetNProtons        (int np) { fNProtons  = np; }
  void SetNNeutrons       (int nn) { fNNeutrons = nn; }
  void SetNSingleGammas   (int ng) { fNSingleGammas = ng ; }
  void SetNRhos           (int nrho_plus, int nrho_0, int nrho_minus);
  void UnsetCharm         (void);
  void UnsetStrange       (void);
  void ResetNPions        (void);
  void ResetNNucleons     (void);
  void ResetNSingleGammas (void) { fNSingleGammas = 0 ;}
  void ResetNRhos         (void);
  void SetResonance       (Resonance_t res);
  void SetDecayMode       (int decay_mode);

  // Copy, reset, print itself and build string code
  void   Reset    (void);                          ///< reset object
  void   Copy     (const XclsTag & xcls);          ///< copy input XclsTag object
  string AsString (void) const;                    ///< pack into a string code
  void   Print    (ostream & stream) const;        ///< print

  XclsTag &        operator =  (const XclsTag & xcls);                  ///< copy
  friend ostream & operator << (ostream& stream, const XclsTag & xcls); ///< print

private:

  // Private data members
  bool        fIsStrangeEvent ;       ///< true if we have strange production
  bool        fIsCharmEvent ;         ///< true if we have charm production
  int         fStrangeHadronPdg;      ///< strange hadron pdg-code
  int         fCharmedHadronPdg;      ///< charmed hadron pdg-code
  int         fNProtons;              ///< # of p's in the hadronic system after this Xcls reaction (before FSI)
  int         fNNeutrons;             ///< # of n's in the hadronic system after this Xcls reaction (before FSI)
  int         fNPi0;                  ///< # of pi^0's in the hadronic system after this Xcls reaction (before FSI)
  int         fNPiPlus;               ///< # of pi^+'s in the hadronic system after this Xcls reaction (before FSI)
  int         fNPiMinus;              ///< # of pi^-'s in the hadronic system after this Xcls reaction (before FSI)
  int         fNSingleGammas;         ///< # of single gammas in the hadronic system after this Xcls reaction (before FSI)
  int         fNRho0;                 ///< # of rho^0's in the hadronic system after this Xcls reaction (before FSI)
  int         fNRhoPlus;              ///< # of rho^+'s in the hadronic system after this Xcls reaction (before FSI)
  int         fNRhoMinus;             ///< # of rho^-'s in the hadronic system after this Xcls reaction (before FSI)
  Resonance_t fResonance;             ///< baryon resonance excited by probe
  int         fDecayMode;

ClassDef(XclsTag,4)
};

}      // namespace

#endif // _FINAL_STATE_H_
