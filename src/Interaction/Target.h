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

using std::ostream;
using std::string;

namespace genie {

class Target {

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

  void SetStruckNucleonPDGCode (int nucl_pdgc);
  void SetStruckNucleonP4      (const TLorentzVector & p4);  

  int              StruckNucleonPDGCode (void) const;
  double           StruckNucleonMass    (void) const;
  TLorentzVector * StruckNucleonP4      (void) const;
  
  double Mass                    (void) const;
  double Charge                  (void) const;
  double BindEnergy              (void) const;
  double BindEnergyPerNucleon    (void) const;
  double BindEnergyLastNucleon   (void) const;
  bool   IsFreeNucleon           (void) const;
  bool   IsProton                (void) const;
  bool   IsNeutron               (void) const;
  bool   IsNucleus               (void) const;
  bool   IsParticle              (void) const;
  bool   IsValidNucleus          (void) const;
  bool   StruckNucleonIsSet      (void) const;
  bool   IsEvenEven              (void) const;
  bool   IsEvenOdd               (void) const;
  bool   IsOddOdd                (void) const;        

//  double        NuclDensityProfile    (double r)              const;
//  NucSpectrum & ExcitationSpectrum    (const)                 const;
//  double        Parity                (void)                  const;
//  double        Lifetime              (void)                  const;
//  double        ElDipoleMoment        (void)                  const;
//  double        ElQuadrupoleMoment    (void)                  const;
//  double        MgDipoleMoment        (void)                  const;
//  double        MeanFreePath          (TParticlePDG & probe)  const;
//  double        Spin                  (void)                  const;
//  double        SpinZ                 (void)                  const;
//  double        Isospin               (void)                  const;
//  double        IsospinZ              (void)                  const;
//  double        OrbitalMomentum       (void)                  const;
//  double        OrbitalMomentumZ      (void)                  const;
//  bool          IsAllowedDecay        (DecayChannel & dc)     const;
//  double        DecayRate             (DecayChannel & dc)     const;
//  DecayList &   GetDecayChannels      (void)                  const;

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
  int  fTgtPDG;
  TLorentzVector * fStruckNucP4;
};

}      // genie namespace

#endif // _TARGET_H_
