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

  void SetZA (int Z, int A);

  //-- Get atomic number, mass number & number of neutrons
  int  Z (void) const;
  int  N (void) const;
  int  A (void) const;

  //-- Get nucleus PDG code according to the MINOS PDG extensions
  int  PDGCode (void) const;

  //-- Set & Get struck nucleon pdg code & 4-momentum

  void SetStruckNucleonPDGCode (int pdgc);
  void SetStruckQuarkPDGCode   (int pdgc);
  void SetStruckNucleonP4      (const TLorentzVector & p4);
  void SetStruckSeaQuark       (bool tf);

  int              StruckNucleonPDGCode (void) const;
  int              StruckQuarkPDGCode   (void) const;
  double           StruckNucleonMass    (void) const;
  TLorentzVector * StruckNucleonP4      (void) const;

  double Mass                    (void) const;
  double Charge                  (void) const;
  bool   IsFreeNucleon           (void) const;
  bool   IsProton                (void) const;
  bool   IsNeutron               (void) const;
  bool   IsNucleus               (void) const;
  bool   IsParticle              (void) const;
  bool   IsValidNucleus          (void) const;
  bool   StruckNucleonIsSet      (void) const;
  bool   StruckQuarkIsSet        (void) const;
  bool   StruckQuarkIsFromSea    (void) const;
  bool   IsEvenEven              (void) const;
  bool   IsEvenOdd               (void) const;
  bool   IsOddOdd                (void) const;

  string AsString (void) const;
  void   Copy     (const Target & t);
  void   Print    (ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const Target & target);

private:

  //-- Initialize
  void Init (void);

  //-- Only valid nucleus & struck nucleon can be set
  void ForceNucleusValidity       (void);
  bool ForceStruckNucleonValidity (void);

  //-- Data members
  int  fZ;
  int  fA;
  int  fStruckNucPDG;
  int  fStruckQuarkPDG;
  int  fTgtPDG;
  bool fStruckSeaQuark;

  TLorentzVector * fStruckNucP4;

ClassDef(Target,1)
};

}      // genie namespace

#endif // _TARGET_H_
