//____________________________________________________________________________
/*!

\program testConfigPool

\brief   test program used for testing / debugging GENIE's AlgConfigPool

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 4, 2004
*/
//____________________________________________________________________________

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/QELFormFactorsModelI.h"
#include "Messenger/Messenger.h"

using namespace genie;

int main(int argc, char ** argv)
{
  //-- Get an instance of the ConfigPool

  LOG("Main", pINFO) << "Get config pool instance";
  AlgConfigPool * pool = ConfigPool::Instance();

  //-- Print the ConfigPool ( => print all its configuration registries )

  LOG("Main", pINFO) << "Printing the config pool";
  LOG("Main", pINFO) << ENDL << *pool;

  //-- instantiate an algorithm

  LOG("Main", pINFO) << "Instantiate a concrete algorithm";
  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * alg_base = algf->GetAlgorithm(
                                  "genie::LlewellynSmithModelCC","Default");
  const QELFormFactorsModelI * llewellyn_smith =
                      dynamic_cast<const QELFormFactorsModelI *> (alg_base);

  //-- Ask the ConfigPool for this algorithm's config registry and print it

  LOG("Main", pINFO) << "Find the configuration registry";

  Registry * config = pool->FindRegistry( llewellyn_smith );

  if(config) LOG("Main", pINFO) << ENDL << *config;

  return 0;
}

