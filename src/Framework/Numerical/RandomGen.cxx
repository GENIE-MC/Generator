//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 24, 2013 - CA
   No longer uses the $GSEED variable for setting the random number seed.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TSystem.h>
#include <TPythia6.h>

#include "Framework/Conventions/Controls.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"

using namespace genie::controls;

namespace genie {

//____________________________________________________________________________
RandomGen * RandomGen::fInstance = 0;
//____________________________________________________________________________
RandomGen::RandomGen()
{
  LOG("Rndm", pINFO) << "RandomGen late initialization";

  fInitalized = false;
  fInstance = 0;
/*
  // try to get this job's random number seed from the environment
  const char * seed = gSystem->Getenv("GSEED");
  if(seed) {
    LOG("Rndm", pDEBUG) << "Reading RandomNumGenerator seed env. var $GSEED";
    fCurrSeed = atoi(seed);
  } else {
    LOG("Rndm", pINFO) << "Env. var. $GSEED is not set. Using default seed";
    fCurrSeed = kDefaultRandSeed; // default seed number
  }
  this->InitRandomGenerators(fCurrSeed);
*/

  if ( gSystem->Getenv("GSEED") ) {
    LOG("Rndm", pFATAL) 
      << "\n\n"
      << "************************************************************************************** \n"
      << "The random number seed is no longer set via the $GSEED variable.\n"
      << "Please use the --seed option implemented in all GENIE apps or, if you access RandomGen \n"
      << "directly in your user code, use RandomGen::SetSeed(long int seed).\n"
      << "Unset $GSEED to continue running GENIE. \n"
      << "************************************************************************************** \n";
    gAbortingInErr = true;
    exit(1);
  }

  fCurrSeed = kDefaultRandSeed; // a default seed number is set a init
  this->InitRandomGenerators(fCurrSeed);

  fInitalized = true;
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
  LOG("Rndm", pNOTICE)
     << "Setting"
     << ((fInitalized) ? " " : " default ")
     << "random number seed" 
     << ((fInitalized) ? ": " : " at random number generator initialization: ")
     << seed;

  // Set the seed number for all internal GENIE random number generators
  this->RndKine ().SetSeed(seed);
  this->RndHadro().SetSeed(seed);
  this->RndDec  ().SetSeed(seed); 
  this->RndFsi  ().SetSeed(seed);
  this->RndLep  ().SetSeed(seed);
  this->RndISel ().SetSeed(seed);
  this->RndGeom ().SetSeed(seed);
  this->RndFlux ().SetSeed(seed);
  this->RndEvg  ().SetSeed(seed);
  this->RndNum  ().SetSeed(seed); 
  this->RndGen  ().SetSeed(seed);

  // Set the seed number for ROOT's gRandom
  gRandom ->SetSeed (seed);

  // Set the PYTHIA6 seed number
  TPythia6 * pythia6 = TPythia6::Instance();
  pythia6->SetMRPY(1, seed);

  LOG("Rndm", pINFO) << "RndKine  seed = " << this->RndKine ().GetSeed();
  LOG("Rndm", pINFO) << "RndHadro seed = " << this->RndHadro().GetSeed();
  LOG("Rndm", pINFO) << "RndDec   seed = " << this->RndDec  ().GetSeed();
  LOG("Rndm", pINFO) << "RndFsi   seed = " << this->RndFsi  ().GetSeed();
  LOG("Rndm", pINFO) << "RndLep   seed = " << this->RndLep  ().GetSeed();
  LOG("Rndm", pINFO) << "RndISel  seed = " << this->RndISel ().GetSeed();
  LOG("Rndm", pINFO) << "RndGeom  seed = " << this->RndGeom ().GetSeed();
  LOG("Rndm", pINFO) << "RndFlux  seed = " << this->RndFlux ().GetSeed();
  LOG("Rndm", pINFO) << "RndEvg   seed = " << this->RndEvg  ().GetSeed();
  LOG("Rndm", pINFO) << "RndNum   seed = " << this->RndNum  ().GetSeed();
  LOG("Rndm", pINFO) << "RndGen   seed = " << this->RndGen  ().GetSeed();
  LOG("Rndm", pINFO) << "gRandom  seed = " << gRandom->GetSeed();
  LOG("Rndm", pINFO) << "PYTHIA6  seed = " << pythia6->GetMRPY(1);
}
//____________________________________________________________________________
void RandomGen::InitRandomGenerators(long int seed)
{
  fRandom3 = new TRandom3();
  this->SetSeed(seed);
}
//____________________________________________________________________________
} // genie namespace
