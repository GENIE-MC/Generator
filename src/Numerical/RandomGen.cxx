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

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TRandom3.h>
#include <TRandom2.h>
#include <TSystem.h>

#include "Conventions/Controls.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"

using namespace genie::controls;

namespace genie {

//____________________________________________________________________________
RandomGen * RandomGen::fInstance = 0;
//____________________________________________________________________________
RandomGen::RandomGen()
{
  LOG("Rndm", pINFO) << "RandomGen late initialization";

  fInstance = 0;

  // try to get this job's random number seed from the environment
  const char * seed = gSystem->Getenv("GSEED");
  if(seed) {
    LOG("Rndm", pDEBUG) << "Reading RandomNumGenerator seed env. var $GSEED";
    fCurrSeed = atoi(seed);
  } else {
    LOG("Rndm", pINFO) << "Env. var. $GSEED is not set. Using default seed";
    fCurrSeed = kDefaultRandSeed; // default seed number
  }
  LOG("Rndm", pINFO)
            << "Starting Random Number Generators with seed = " << fCurrSeed;

  this->InitRandomGenerators(fCurrSeed);
}
//____________________________________________________________________________
RandomGen::~RandomGen()
{
  fInstance = 0;
}
//____________________________________________________________________________
RandomGen * RandomGen::Instance()
{
  if(fInstance == 0) {
    static RandomGen::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new RandomGen;
  }
  return fInstance;
}
//____________________________________________________________________________
void RandomGen::SetSeed(long int seed)
{
  gRandom  -> SetSeed (seed);
  fRandom1 -> SetSeed (seed);
  fRandom2 -> SetSeed (seed);
  fRandom3 -> SetSeed (seed);
}
//____________________________________________________________________________
void RandomGen::InitRandomGenerators(long int seed)
{
  gRandom->SetSeed(seed);

  fRandom1 = new TRandom  (seed);
  fRandom2 = new TRandom2 (seed);
  fRandom3 = new TRandom3 (seed);
}
//____________________________________________________________________________
} // genie namespace
