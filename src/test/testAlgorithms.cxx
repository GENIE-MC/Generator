//____________________________________________________________________________
/*!

\program testAlgorithms

\brief   program used for testing / debugging GENIE's Algorithms & AlgFactory

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 26, 2006

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/QELFormFactorsModelI.h"
#include "Base/ELFormFactorsModelI.h"
#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"

using namespace genie;

void testReconfigInCommonPool   (void);
void testReconfigInOwnedModules (void);

int main(int /*argc*/, char ** /*argv*/)
{
  testReconfigInCommonPool();
  testReconfigInOwnedModules();

  return 0;
}
//____________________________________________________________________________
void testReconfigInCommonPool(void)
{
// Test reconfiguration of algorithms in common pool 
// (AP: algorithm factory pool, CP: configuration pool)
//
// The test function access an algorithm stored in the AP & then access the 
// configuration registry that the algorithm is looking up at the CP. 
// Then it modifies the configuration stored at the CP, and forces all 
// algorithms stored at the AP to reconfigure themselves.
//
// Reconfiguring common pool objects makes it easier to guarantee consistency:
// Eg say that you have a low-level algorithm stored at the AP that is 
// referenced by many higher-level objects, also stored at the AP.
// Reconfiguring this instance of the low-level algorithms **automatically**
// guarantees that all high-level algorithms use the reconfigured version.
// This does *not* require you to know which high-level algorithms are
// pointing to which low-level algorithms which is typically difficult to
// keep track of

  // get the algorithm factory & config pool
  AlgConfigPool * cnfp = AlgConfigPool::Instance();
  AlgFactory *    algf = AlgFactory::Instance();

  // instantiate some algorithms
  LOG("Main", pINFO) << "Instantiate a concrete algorithm";
  AlgId id("genie::QELPXSec","CC-Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  LOG("Main", pINFO) << *alg;

  LOG("Main", pINFO) << "Access its configuration at the config pool";
  Registry * r = cnfp->FindRegistry(alg);
  r->UnLock();
  LOG("Main", pINFO) << *r;

  // modify configuration 
  LOG("Main", pINFO) << "Modifying registry";
  r->Set("CabbiboAngle",0.25);
  LOG("Main", pINFO) << *r;

  // force reconfiguration
  algf->ForceReconfiguration();

  // print all algorithms stored at the algorithm factory pool (for now just 1) 
  // & their configurations to verify that the change was propagated correctly
  LOG("Main", pINFO) << *algf;
}
//____________________________________________________________________________
void testReconfigInOwnedModules (void)
{
// Test reconfiguration of owned "algorithm/configuration blocks"
//
// In GENIE it is possible to take ownership of an algorithms and its 
// subtructure (all the sub-algorithms it is referencing) by extracting (clones 
// of) them from the AlgFactory pool to a local pool. The owned algorithms are 
// also forced to take ownership of their configurations. All the configuration 
// variables are bundled together at the top level algorithms' configuration.
//
// Having a group of algorithms/configurations behaving as a monolithic block, 
// with a single point of configuration (the top level) is to be used when bits 
// & pieces of GENIE are used in isolation for data fitting or reweighting
//
// For more details and information on the naming convention used when all
// sub-algorithms configurations are merged to the top level algorithm config
// see the Algorithm package
//

  // get the algorithm factory & config pool
  AlgFactory * algf = AlgFactory::Instance();

  // instantiate some algorithms
  LOG("Main", pINFO) << "Instantiate a concrete algorithm";
  AlgId id("genie::QELPXSec","CC-Default");
  Algorithm * alg = algf->AdoptAlgorithm(id);
  XSecAlgorithmI * xsecalg = dynamic_cast<XSecAlgorithmI*>(alg);

  LOG("Main", pINFO) << *xsecalg;

  LOG("Main", pINFO) << "Adopting substructure";
  xsecalg->AdoptSubstructure();

  LOG("Main", pINFO) << *xsecalg;

  // access the top level algorithm registry where all the config params
  // (including config params of all the referenced sub-algs have been 
  // bundled) following a special naming convention
  //
  LOG("Main", pINFO) << "Taking a clone of the deep config registry:";
  Registry r(xsecalg->GetConfig());
  LOG("Main", pINFO) << r;

  LOG("Main", pINFO) << "Modifying parameters at the deep registry";
  r.Set("CabbiboAngle",                        0.23); // refers to top level alg
  r.Set("FormFactorsAlg/MuN",                 -1.92); // refers to an alg 1-level deep
  r.Set("FormFactorsAlg/ElFormFactorsAlg/MuN",-1.92); // refers to an alg 2-level deep

  LOG("Main", pINFO) << "Modified deep config registry:";
  LOG("Main", pINFO) << r;

  // This would reconfigure the top level algorithm and then, in a recursive
  // mode, all owned sub-algorithms will be reconfigured and then their owned
  // sub-algorithms and so on, however complex the algorithm strucure
  xsecalg->Configure(r);

  LOG("Main", pINFO) << *xsecalg;
}
//____________________________________________________________________________
