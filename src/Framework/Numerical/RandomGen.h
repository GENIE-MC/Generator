//____________________________________________________________________________
/*!

\class    genie::RandomGen

\brief    A singleton holding random number generator classes. All random 
          number generation in GENIE should take place through this class.
          Ensures that the random number generator seed is set consistently
          to all GENIE modules and that all modules use the preferred rndm
          number generator.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  September 22, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RANDOM_GEN_H_
#define _RANDOM_GEN_H_

#include <TRandom3.h>

namespace genie {

class RandomGen {

public:

  //! Access instance
  static RandomGen * Instance();

  //! Random number generators used by various GENIE modules.
  //! (See note at http://root.cern.ch/root/html//TRandom.html
  //!  on using several TRandom objects each with each own
  //!  "independent" run sequence).

  //! At this point, since the actual random number generator
  //! periodicity is very high, all the generators are in fact one!
  //! However, the option to use many generators is reserved.

  //! Currently, the preferred generator is the "Mersenne Twister"
  //! with a periodicity of 10**6000
  //! See: http://root.cern.ch/root/html/TRandom3.html

  //! rnd number generator used by kinematics generators
  TRandom3 & RndKine (void) const { return *fRandom3; } 

  //! rnd number generator used by hadronization models 
  TRandom3 & RndHadro (void) const { return *fRandom3; }

  //! rnd number generator used by decay models 
  TRandom3 & RndDec (void) const { return *fRandom3; }

  //! rnd number generator used by intranuclear cascade monte carlos
  TRandom3 & RndFsi (void) const { return *fRandom3; }

  //! rnd number generator used by final state primary lepton generators
  TRandom3 & RndLep (void) const { return *fRandom3; } 

  //! rnd number generator used by interaction selectors
  TRandom3 & RndISel (void) const { return *fRandom3; }

  //! rnd number generator used by geometry drivers
  TRandom3 & RndGeom (void) const { return *fRandom3; }

  //! rnd number generator used by flux drivers
  TRandom3 & RndFlux (void) const { return *fRandom3; }

  //! rnd number generator used by the event generation drivers
  TRandom3 & RndEvg (void) const { return *fRandom3; }

  //! rnd number generator used by MC integrators & other numerical methods
  TRandom3 & RndNum (void) const { return *fRandom3; }

  //! rnd number generator for generic usage
  TRandom3 & RndGen  (void) const { return *fRandom3; }

  long int GetSeed (void)         const { return fCurrSeed; }
  void     SetSeed (long int seed);

private:

  RandomGen();
  RandomGen(const RandomGen & rgen);
  virtual ~RandomGen();

  static RandomGen * fInstance;

  TRandom3 * fRandom3;    ///< Mersenne Twistor
  long int   fCurrSeed;   ///< random number generator seed number
  bool       fInitalized; ///< done initializing singleton?

  void InitRandomGenerators(long int seed);

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (RandomGen::fInstance !=0) {
            delete RandomGen::fInstance;
            RandomGen::fInstance = 0;
         }
      }
  };

  friend struct Cleaner;
};

}      // genie namespace

#endif // _RANDOM_GEN_H_
