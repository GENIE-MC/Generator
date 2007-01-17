//____________________________________________________________________________
/*!

\program testConfigPool

\brief   test program used for testing / debugging GENIE's AlgConfigPool

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/QELFormFactorsModelI.h"
#include "Base/ELFormFactorsModelI.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  //-- Get an instance of the ConfigPool

  LOG("Main", pINFO) << "Get config pool instance";
  AlgConfigPool * pool = AlgConfigPool::Instance();

  //-- Print the ConfigPool ( => print all its configuration registries )

  LOG("Main", pINFO) << "Printing the config pool";
  LOG("Main", pINFO) << ENDL << *pool;

  //-- get the algorithm factory

  AlgFactory * algf = AlgFactory::Instance();

  //-- instantiate an algorithm

  LOG("Main", pINFO) << "Instantiate a concrete algorithm";

  const Algorithm * alg0 = algf->GetAlgorithm(
                                  "genie::LlewellynSmithModelCC","Default");
  const QELFormFactorsModelI * llewellyn_smith =
                         dynamic_cast<const QELFormFactorsModelI *> (alg0);
  LOG("Main", pINFO) << *alg0;


  LOG("Main", pINFO) << "Instantiate another concrete algorithm";

  const Algorithm * alg1 = algf->GetAlgorithm(
                             "genie::DipoleELFormFactorsModel","Default");
  const ELFormFactorsModelI * dipole_elff =
                         dynamic_cast<const ELFormFactorsModelI *> (alg1);
  LOG("Main", pINFO) << *alg1;


  //-- Ask the ConfigPool for this algorithm's config registry and print it

  LOG("Main", pINFO) << "Find the configuration for both algorithms";

  Registry * config1 = pool->FindRegistry( llewellyn_smith );
  Registry * config2 = pool->FindRegistry( dipole_elff     );

  if(config1) LOG("Main", pINFO) << ENDL << *config1;
  if(config2) LOG("Main", pINFO) << ENDL << *config2;

  return 0;
}

