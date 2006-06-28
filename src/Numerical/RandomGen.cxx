//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - September 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

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
  if(fRandom3) delete fRandom3;
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
  fRandom3 -> SetSeed (seed);

  LOG("Rndm", pINFO) << "gRandom  seed = " << gRandom  -> GetSeed();
  LOG("Rndm", pINFO) << "fRandom3 seed = " << fRandom3 -> GetSeed();
}
//____________________________________________________________________________
void RandomGen::InitRandomGenerators(long int seed)
{
  fRandom3 = new TRandom3();

  this->SetSeed(seed);
}
//____________________________________________________________________________
} // genie namespace
