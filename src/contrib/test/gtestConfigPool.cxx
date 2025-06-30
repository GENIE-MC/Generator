//____________________________________________________________________________
/*!

\program gtestConfigPool

\brief   Program used for testing / debugging GENIE's AlgConfigPool

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created May 4, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  // Get an instance of the ConfigPool
  LOG("test", pINFO) << "Get config pool instance";
  AlgConfigPool * pool = AlgConfigPool::Instance();

  // Print the ConfigPool ( => print all its configuration registries )
  LOG("test", pINFO) << "Printing the config pool\n" << *pool;

  // Get the algorithm factory
  AlgFactory * algf = AlgFactory::Instance();

  // Instantiate an algorithm
  LOG("test", pINFO) << "Instantiate a concrete algorithm";
  const Algorithm * alg0 = 
       algf->GetAlgorithm("genie::LwlynSmithFFCC","Default");
  const QELFormFactorsModelI * llewellyn_smith =
       dynamic_cast<const QELFormFactorsModelI *> (alg0);
  LOG("test", pINFO) << *alg0;

  LOG("test", pINFO) << "Instantiate another concrete algorithm";
  const Algorithm * alg1 = 
       algf->GetAlgorithm("genie::DipoleELFormFactorsModel","Default");
  const ELFormFactorsModelI * dipole_elff =
        dynamic_cast<const ELFormFactorsModelI *> (alg1);
  LOG("test", pINFO) << *alg1;

  // Ask the ConfigPool for this algorithm's config registry and print it

  LOG("test", pINFO) << "Find the configuration for both algorithms";

  Registry * config1 = pool->FindRegistry( llewellyn_smith );
  Registry * config2 = pool->FindRegistry( dipole_elff     );

  if(config1) LOG("test", pINFO) << "1st algorithm config: \n" << *config1;
  if(config2) LOG("test", pINFO) << "2nd algorithm config: \n" << *config2;

  return 0;
}

