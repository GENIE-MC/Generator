//____________________________________________________________________________
/*!

\class    genie::Target

\brief    A Neutrino Interaction Target. Is a transparent encapsulation of
          quite different physical systems such as a nuclear target, a
          'spectator' nuclear target with a Hit nucleon, a free nucleon or
          a free particle (eg a e- target in the inverse muon decay reaction)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 03, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
  Target(int Z, int A, int hit_nucleon_pdgc);
  Target(const Target & tgt);
  Target(TRootIOCtor*);
 ~Target();

  //-- Set target properties

  void SetId                  (int pdgc);
  void SetId                  (int Z, int A);
  void SetHitNucPdg           (int pdgc);
  void SetHitNucP4            (const TLorentzVector & p4);
  void SetHitNucPosition        (double r);
  void SetHitQrkPdg           (int pdgc);
  void SetHitSeaQrk           (bool tf);
  void ForceHitNucOnMassShell (void);

  //-- Query target information

  int    Z              (void) const { return fZ;      }
  int    N              (void) const { return fA-fZ;   }
  int    A              (void) const { return fA;      }
  int    Pdg            (void) const { return fTgtPDG; }
  double Mass           (void) const;
  double Charge         (void) const;
  bool   IsFreeNucleon  (void) const;
  bool   IsProton       (void) const;
  bool   IsNeutron      (void) const;
  bool   IsNucleus      (void) const;
  bool   IsParticle     (void) const;
  bool   IsValidNucleus (void) const;
  bool   HitNucIsSet    (void) const;
  bool   HitQrkIsSet    (void) const;
  bool   HitSeaQrk      (void) const;
  bool   IsEvenEven     (void) const;
  bool   IsEvenOdd      (void) const;
  bool   IsOddOdd       (void) const;
  int    HitNucPdg      (void) const;
  int    HitQrkPdg      (void) const;
  double HitNucMass     (void) const;
  double HitNucPosition (void) const { return fHitNucRad; }

  const TLorentzVector & HitNucP4    (void) const { return *this->HitNucP4Ptr(); }
  TLorentzVector *       HitNucP4Ptr (void) const;
  
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
  void ForceNucleusValidity (void);
  bool ForceHitNucValidity  (void);
  void AutoSetHitNuc        (void);

  //-- Private data members
  int  fZ;                    ///< nuclear target Z
  int  fA;                    ///< nuclear target A
  int  fTgtPDG;               ///< nuclear target PDG code
  int  fHitNucPDG;            ///< hit nucleon PDG code
  int  fHitQrkPDG;            ///< hit quark PDG code
  bool fHitSeaQrk;            ///< hit quark from sea?
  TLorentzVector * fHitNucP4; ///< hit nucleon 4p
  double fHitNucRad;          ///< hit nucleon position

ClassDef(Target,2)
};

}      // genie namespace

#endif // _TARGET_H_

