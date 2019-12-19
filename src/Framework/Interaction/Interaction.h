//____________________________________________________________________________
/*!

\class    genie::Interaction

\brief    Summary information for an interaction.

          It is a container of an InitialState, a ProcessInfo, an XclsTag
          and a Kinematics object.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Changes required to implement the GENIE Boosted Dark Matter module 
          were installed by Josh Berger (Univ. of Wisconsin)

\created  April 25, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <ostream>
#include <string>

#include <TObject.h>

#include "Framework/Conventions/RefFrame.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/Kinematics.h"
#include "Framework/Interaction/XclsTag.h"
#include "Framework/Interaction/KPhaseSpace.h"

using std::ostream;
using std::string;

class TRootIOCtor;

namespace genie {

const UInt_t kISkipProcessChk      = 1<<17; ///< if set, skip process validity checks
const UInt_t kISkipKinematicChk    = 1<<16; ///< if set, skip kinematic validity checks
const UInt_t kIAssumeFreeNucleon   = 1<<15; ///<
const UInt_t kIAssumeFreeElectron  = 1<<15; ///<
const UInt_t kINoNuclearCorrection = 1<<14; ///< if set, inhibit nuclear corrections 

class Interaction;
ostream & operator << (ostream & stream, const Interaction & i); 

class Interaction : public TObject {

public:
  using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TObject::Copy;  // 

  Interaction();
  Interaction(const InitialState & init, const ProcessInfo & proc);
  Interaction(const Interaction & i);
  Interaction(TRootIOCtor*);
 ~Interaction();

  // Methods accessing aggregate/owned objects holding interaction information
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

  // Methods to set interaction's properties
  void SetInitState (const InitialState & init);
  void SetProcInfo  (const ProcessInfo &  proc);
  void SetKine      (const Kinematics &   kine);
  void SetExclTag   (const XclsTag &      xcls);

  // Get the final state primary lepton and recoil nucleon (if) uniquely
  // determined for the specified interaction
  int            FSPrimLeptonPdg  (void) const; ///< final state primary lepton pdg
  int            RecoilNucleonPdg (void) const; ///< recoil nucleon pdg
  TParticlePDG * FSPrimLepton     (void) const; ///< final state primary lepton
  TParticlePDG * RecoilNucleon    (void) const; ///< recoil nucleon 

  // Copy, reset, print itself and build string code
  void   Reset    (void);
  void   Copy     (const Interaction & i);
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  // Overloaded operators
  Interaction &    operator =  (const Interaction & i);                   ///< copy
  friend ostream & operator << (ostream & stream, const Interaction & i); ///< print

  // Use the "Named Constructor" C++ idiom for fast creation of typical interactions
  static Interaction * DISCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DISCC     (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DISCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DISCC     (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);
  static Interaction * DISNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DISNC     (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DISNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DISNC     (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);
  static Interaction * DISEM     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DISEM     (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DISEM     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DISEM     (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);
  static Interaction * QELCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * QELCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * QELNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * QELNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * QELEM     (int tgt, int nuc, int probe, double E=0);
  static Interaction * QELEM     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * IBD       (int tgt, int nuc, int probe, double E=0);
  static Interaction * IBD       (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * RESCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * RESCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * RESNC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * RESNC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * RESEM     (int tgt, int nuc, int probe, double E=0);
  static Interaction * RESEM     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DFRCC     (int tgt, int nuc, int probe, double E=0);
  static Interaction * DFRCC     (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * COHCC     (int tgt, int probe, unsigned int prod_pdg, double E=0);
  static Interaction * COHCC     (int tgt, int probe, unsigned int prod_pdg, 
				  const TLorentzVector & p4probe);
  static Interaction * COHNC     (int tgt, int probe, unsigned int prod_pdg, double E=0);
  static Interaction * COHNC     (int tgt, int probe, unsigned int prod_pdg, 
				  const TLorentzVector & p4probe);
  static Interaction * CEvNS     (int tgt, int probe, double E=0);
  static Interaction * CEvNS     (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * IMD       (int tgt, double E=0);
  static Interaction * IMD       (int tgt, const TLorentzVector & p4probe);
  static Interaction * AMNuGamma (int tgt, int nuc, int probe, double E=0);
  static Interaction * AMNuGamma (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * MECCC     (int tgt, int nuccluster, int probe, double E=0);
  static Interaction * MECCC     (int tgt, int nuccluster, int probe, const TLorentzVector & p4probe);
  static Interaction * MECCC     (int tgt, int probe, double E=0);
  static Interaction * MECCC     (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * MECNC     (int tgt, int nuccluster, int probe, double E=0);
  static Interaction * MECNC     (int tgt, int nuccluster, int probe, const TLorentzVector & p4probe);
  static Interaction * MECEM     (int tgt, int nuccluster, int probe, double E=0);
  static Interaction * MECEM     (int tgt, int nuccluster, int probe, const TLorentzVector & p4probe);
  static Interaction * GLR       (int tgt, double E=0);
  static Interaction * GLR       (int tgt, const TLorentzVector & p4probe);
  static Interaction * NDecay    (int tgt, int decay_mode=-1, int decayed_nucleon = 0);
  static Interaction * NOsc      (int tgt, int annihilation_mode=-1);
  static Interaction * ASK       (int tgt, int probe, double E=0);
  static Interaction * ASK       (int tgt, int probe, const TLorentzVector & p4probe);
  static Interaction * DME       (int tgt, int nuc, int probe, double E=0);
  static Interaction * DME       (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DMDI      (int tgt, int nuc, int probe, double E=0);
  static Interaction * DMDI      (int tgt, int nuc, int qrk, bool sea, int probe, double E=0);
  static Interaction * DMDI      (int tgt, int nuc, int probe, const TLorentzVector & p4probe);
  static Interaction * DMDI      (int tgt, int nuc, int qrk, bool sea, int probe, const TLorentzVector & p4probe);

private:

  // Methods for Interaction initialization and clean up
  void Init    (void);
  void CleanUp (void);

  // Utility method for "named ctor"
  static Interaction * Create(int tgt, int probe, ScatteringType_t st, InteractionType_t it);

  // Private data members
  InitialState * fInitialState;  ///< Initial State info
  ProcessInfo *  fProcInfo;      ///< Process info (scattering, weak current,...)
  Kinematics *   fKinematics;    ///< kinematical variables
  XclsTag *      fExclusiveTag;  ///< Additional info for exclusive channels
  KPhaseSpace *  fKinePhSp;      ///< Kinematic phase space
  
ClassDef(Interaction,2)
};

}      // genie namespace

#endif // _INTERACTION_H_
