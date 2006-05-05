//____________________________________________________________________________
/*!

\class    genie::Target

\brief    A Neutrino Interaction Target. Is a transparent encapsulation of
          quite different physical systems such as a nuclear target, a
          'spectator' nuclear target with a struck nucleon, a free nucleon or
          a free particle (eg a e- target in the inverse muon decay reaction)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
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

namespace genie {

class Target : public TObject {

public:
  Target();
  Target(int pdgc);
  Target(int Z, int A);
  Target(int Z, int A, int struck_nucleon_pdgc);
  Target(const Target & tgt);
  ~Target();

  //! Set target properties
  void SetId                         (int pdgc);
  void SetId                         (int Z, int A);
  void SetStruckNucleonPDGCode       (int pdgc);
  void SetStruckNucleonP4            (const TLorentzVector & p4);
  void SetStruckQuarkPDGCode         (int pdgc);
  void SetStruckSeaQuark             (bool tf);
  void ForceStruckNucleonOnMassShell (void);

  //! Get atomic number, mass number & number of neutrons,
  //! nucleus PDG code according to the MINOS PDG extensions
  //! and provide shortcuts for getting the mass/charge
  int    Z       (void) const { return fZ;      }
  int    N       (void) const { return fA-fZ;   }
  int    A       (void) const { return fA;      }
  int    PDGCode (void) const { return fTgtPDG; }

  //! Query for target information
  double Mass                 (void) const;
  double Charge               (void) const;
  bool   IsFreeNucleon        (void) const;
  bool   IsProton             (void) const;
  bool   IsNeutron            (void) const;
  bool   IsNucleus            (void) const;
  bool   IsParticle           (void) const;
  bool   IsValidNucleus       (void) const;
  bool   StruckNucleonIsSet   (void) const;
  bool   StruckQuarkIsSet     (void) const;
  bool   StruckQuarkIsFromSea (void) const;
  bool   IsEvenEven           (void) const;
  bool   IsEvenOdd            (void) const;
  bool   IsOddOdd             (void) const;
  int    StruckNucleonPDGCode (void) const;
  int    StruckQuarkPDGCode   (void) const;
  double StruckNucleonMass    (void) const;
  TLorentzVector * StruckNucleonP4 (void) const;

  //! Copy, reset, compare, print itself and build string code
  void   Reset    (void);
  void   Copy     (const Target & t);
  bool   Compare  (const Target & t) const;
  string AsString (void) const;
  void   Print    (ostream & stream) const;

  bool             operator == (const Target & t) const;
  Target &         operator =  (const Target & t);
  friend ostream & operator << (ostream & stream, const Target & t);

private:

  //! Methods for Target initialization and clean up
  void Init    (void);
  void CleanUp (void);

  //! Methods assuring nucleus & struck nucleon validity
  void ForceNucleusValidity       (void);
  bool ForceStruckNucleonValidity (void);
  void AutoSetStruckNucleon       (void);

  //! Private data members
  int  fZ;
  int  fA;
  int  fTgtPDG;
  int  fStruckNucPDG;
  int  fStruckQuarkPDG;
  bool fStruckSeaQuark;
  TLorentzVector * fStruckNucP4;

ClassDef(Target,1)
};

}      // genie namespace

#endif // _TARGET_H_

