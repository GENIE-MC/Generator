//____________________________________________________________________________
/*!

\class    genie::RandomGen

\brief    A singleton holding ROOT's random number generator classes
          (TRandom, TRandom2, TRandom3).

          All random number generation in GENIE should take place through
          this class. \n

          It prevents random number generator objects from being created (and
          therefore re-initialized) by the different GENIE event generation
          objects. In that case, one could end-up using the same few numbers
          of a random number sequence again and again... \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 22, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _RANDOM_GEN_H_
#define _RANDOM_GEN_H_

#include <TRandom.h>
#include <TRandom2.h>
#include <TRandom3.h>

namespace genie {

class RandomGen {

public:

  static RandomGen * Instance();

  // Simple random number generator
  TRandom  & Random1 (void) const { return *fRandom1; }

  // Generator with periodicity > 10e14
  TRandom2 & Random2 (void) const { return *fRandom2; }

  // Mersenne Twistor
  TRandom3 & Random3 (void) const { return *fRandom3; }

  long int GetSeed (void)         const { return fCurrSeed; }
  void     SetSeed (long int seed);

private:

  RandomGen();
  RandomGen(const RandomGen & rgen);
  virtual ~RandomGen();

  static RandomGen * fInstance;

  TRandom  * fRandom1; // Simple random number generator
  TRandom2 * fRandom2; // Generator with periodicity > 10e14
  TRandom3 * fRandom3; // Mersenne Twistor

  long int fCurrSeed;

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
