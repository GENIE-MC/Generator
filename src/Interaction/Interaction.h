//____________________________________________________________________________
/*!

\class    genie::Interaction

\brief    Summary information for an interaction.

          It is a container of an InitialState, a ProcessInfo, an XclsTag
          and a Kinematics object.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  April 25, 2004

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <ostream>
#include <string>

#include <TObject.h>

#include "Conventions/RefFrame.h"
#include "Interaction/InitialState.h"
#include "Interaction/ProcessInfo.h"
#include "Interaction/Kinematics.h"
#include "Interaction/XclsTag.h"
#include "Interaction/KPhaseSpace.h"

using std::ostream;
using std::string;

namespace genie {

const UInt_t kISkipProcessChk      = 1<<17;
const UInt_t kISkipKinematicChk    = 1<<16;
const UInt_t kIAssumeFreeNucleon   = 1<<15;
const UInt_t kIAssumeFreeElectron  = 1<<15;
const UInt_t kINoNuclearCorrection = 1<<14;

class Interaction : public TObject {

public:
  Interaction();
  Interaction(const InitialState & init, const ProcessInfo & proc);
  Interaction(const Interaction & i);
 ~Interaction();

  //-- Methods accessing aggregate/owned objects holding interaction information
  const InitialState & InitState     (void) const { return *fInitialState; }
  const ProcessInfo &  ProcInfo      (void) const { return *fProcInfo;     }
  const Kinematics &   Kine          (void) const { return *fKinematics;   }
  const XclsTag &      ExclTag       (void) const { return *fExclusiveTag; }
  const KPhaseSpace &  PhaseSpace    (void) const { return *fKinePhSp;     }
  InitialState *       InitStatePtr  (void) const { return fInitialState;  }
  ProcessInfo *        ProcInfoPtr   (void) const { return fProcInfo;      }
  Kinematics *         KinePtr       (void) const { return fKinematics;    }
  XclsTag *            ExclTagPtr    (void) const { return fExclusiveTag;  }
  KPhaseSpace *        PhaseSpacePtr (void) const { return fKinePhSp;      }

  //-- Methods to set interaction's properties
  void SetInitState (const InitialState & init);
  void SetProcInfo  (const ProcessInfo &  proc);
  void SetKine      (const Kinematics &   kine);
  void SetExclTag   (const XclsTag &      xcls);

  //-- Get the final state primary lepton and recoil nucleon (if) uniquely
  //   determined for the specified interaction
  int            FSPrimLeptonPdg  (void) const; ///< final state primary lepton pdg
  int            RecoilNucleonPdg (void) const; ///< recoil nucleon pdg
  TParticlePDG * FSPrimLepton     (void) const; ///< final state primaru lepton
  TParticlePDG * RecoilNucleon    (void) const; ///< recoil nucleon 

  //-- Copy, reset, print itself and build string code
  void   Reset    (void);
  void   Copy     (const Interaction & i);
  bool   Compare  (const Interaction & i) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  //-- Overloaded operators
  bool             operator == (const Interaction & i) const;             ///< compare
  Interaction &    operator =  (const Interaction & i);                   ///< copy
  friend ostream & operator << (ostream & stream, const Interaction & i); ///< print

  //-- Use the "Named Constructor" C++ idiom for fast creation of typical interactions
  static Interaction * DISCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DISCC     (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DISCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DISCC     (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);
  static Interaction * DISNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DISNC     (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DISNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DISNC     (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);
  static Interaction * QELCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * QELCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * QELNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * QELNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * RESCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * RESCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * RESNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * RESNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DFRCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DFRCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * COHPiCC   (int tgt, int probe, double E=0);
  static Interaction * COHPiCC   (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * COHPiNC   (int tgt, int probe, double E=0);
  static Interaction * COHPiNC   (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * COHEl     (int tgt, int probe, double E=0);
  static Interaction * COHEl     (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * IMD       (int tgt, double E=0);
  static Interaction * IMD       (int tgt, const TLorentzVector & p4probe);
  static Interaction * AMNuGamma (int tgt, int nuc, int probe, double E=0);
  static Interaction * AMNuGamma (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * MECCC     (int tgt, int probe, double E=0);
  static Interaction * MECCC     (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * MECNC     (int tgt, int probe, double E=0);
  static Interaction * MECNC     (int tgt, int probe, const TLorentzVector & p4probe);

private:

  //-- Methods for Interaction initialization and clean up
  void Init    (void);
  void CleanUp (void);

  //-- Utility method for "named ctor"
  static Interaction * Create(int tgt, int probe, ScatteringType_t st, InteractionType_t it);

  //-- Private data members
  InitialState * fInitialState;  ///< Initial State info
  ProcessInfo *  fProcInfo;      ///< Process info (scattering, weak current,...)
  Kinematics *   fKinematics;    ///< kinematical variables
  XclsTag *      fExclusiveTag;  ///< Additional info for exclusive channels
  KPhaseSpace *  fKinePhSp;      ///< Kinematic phase space
  
ClassDef(Interaction,1)
};

}      // genie namespace

#endif // _INTERACTION_H_
