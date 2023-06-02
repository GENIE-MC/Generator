//____________________________________________________________________________
/*!

\class    genie::Target

\brief    A Neutrino Interaction Target. Is a transparent encapsulation of
          quite different physical systems such as a nuclear target, a
          'spectator' nuclear target with a Hit nucleon, a free nucleon or
          a free particle (eg a e- target in the inverse muon decay reaction)

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

          Changes required to implement the Electron Velocity module
          were installed by Brinden Carlson (Univ. of Florida)

\created  May 03, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _TARGET_H_
#define _TARGET_H_

#include <ostream>
#include <string>

#include <TLorentzVector.h>
#include <TObject.h>

using std::ostream;
using std::string;

class TRootIOCtor;

namespace genie {

class Target;
ostream & operator << (ostream & stream, const Target & t);

class Target : public TObject {

using TObject::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
using TObject::Compare;
using TObject::Copy;

public:
  Target();
  Target(int pdgc);
  Target(int Z, int A);
  Target(int Z, int A, int hit_particle_pdgc);
  Target(const Target & tgt);
  Target(TRootIOCtor*);
 ~Target();

  //-- Set target properties

<<<<<<< HEAD
  void SetId                  (int pdgc);
  void SetId                  (int Z, int A);
  void SetHitNucPdg           (int pdgc);
  void SetHitNucP4            (const TLorentzVector & p4);
  void SetHitEleP4            (const TLorentzVector & p4);
  void SetHitNucPosition        (double r);
  void SetHitQrkPdg           (int pdgc);
  void SetHitSeaQrk           (bool tf);
  void ForceHitNucOnMassShell (void);
=======
  void SetId                   (int pdgc);
  void SetId                   (int Z, int A);
  void SetHitPartPdg           (int pdgc);
  void SetHitPartP4            (const TLorentzVector & p4);
  void SetHitPartPosition      (double r);
  void SetHitQrkPdg            (int pdgc);
  void SetHitSeaQrk            (bool tf);
  void ForceHitPartOnMassShell (void);
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb

  //-- Query target information

  int    Z              (void) const { return fZ;      }
  int    N              (void) const { return fA-fZ;   }
  int    A              (void) const { return fA;      }
  int    Pdg            (void) const { return fTgtPDG; }
  double Mass           (void) const;
  double Charge         (void) const;
  bool   IsFreeNucleon  (void) const;
  bool   IsProton       (void) const;
  bool   IsElectron     (void) const;
  bool   IsNeutron      (void) const;
  bool   IsNucleus      (void) const;
  bool   IsParticle     (void) const;
  bool   IsValidNucleus (void) const;
<<<<<<< HEAD
  bool   HitNucIsSet    (void) const;
  bool   HitEleIsSet    (void) const;
=======
  bool   HitPartIsSet   (void) const;
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb
  bool   HitQrkIsSet    (void) const;
  bool   HitSeaQrk      (void) const;
  bool   IsEvenEven     (void) const;
  bool   IsEvenOdd      (void) const;
  bool   IsOddOdd       (void) const;
  int    HitPartPdg     (void) const;
  int    HitQrkPdg      (void) const;
  double HitPartMass    (void) const;
  double HitPartPosition(void) const { return fHitPartRad; }

<<<<<<< HEAD
  const TLorentzVector & HitNucP4    (void) const { return *this->HitNucP4Ptr(); }
  TLorentzVector *       HitNucP4Ptr (void) const;
  const TLorentzVector & HitEleP4    (void) const { return *this->HitEleP4Ptr(); }
  TLorentzVector *       HitEleP4Ptr (void) const;
=======
  const TLorentzVector & HitPartP4    (void) const { return *this->HitPartP4Ptr(); }
  TLorentzVector *       HitPartP4Ptr (void) const;
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb

  //-- Copy, reset, compare, print itself and build string code
  void   Reset    (void);
  void   Copy     (const Target & t);
  bool   Compare  (const Target & t) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const Target & t) const;             ///< equal?
  Target &         operator =  (const Target & t);                   ///< copy
  friend ostream & operator << (ostream & stream, const Target & t); ///< print

private:

  //-- Methods for Target initialization and clean up
  void Init    (void);
  void CleanUp (void);

  //-- Methods assuring nucleus & hit nucleon validity
  void ForceNucleusValidity  (void);
  bool ForceHitPartValidity  (void);
  void AutoSetHitPart        (void);

  //-- Private data members
  int  fZ;                    ///< nuclear target Z
  int  fA;                    ///< nuclear target A
  int  fTgtPDG;               ///< nuclear target PDG code
  int  fHitPartPDG;           ///< hit particle PDG code
  int  fHitQrkPDG;            ///< hit quark PDG code
  bool fHitSeaQrk;            ///< hit quark from sea?
<<<<<<< HEAD
  TLorentzVector * fHitNucP4; ///< hit nucleon 4p
  TLorentzVector * fHitEleP4; ///< hit electron 4p - changes index
  double fHitNucRad;          ///< hit nucleon position
=======
  TLorentzVector *fHitPartP4; ///< hit particle 4p
  double fHitPartRad;         ///< hit particle position
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb

ClassDef(Target,3)
};

}      // genie namespace

#endif // _TARGET_H_
