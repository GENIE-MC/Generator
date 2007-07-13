//____________________________________________________________________________
/*!

\program testDecay

\brief   test program used for testing / debugging the PYTHIA/JETSET
         interface and the GENIE DecayModelI algorithms

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <ostream>

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TClonesArray.h>
#include <TIterator.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Decay/DecayModelI.h"
#include "Decay/PythiaDecayer.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;

using std::ostream;
using std::endl;

ostream & operator<< (ostream & stream, const TClonesArray * particle_list);
ostream & operator<< (ostream & stream, const TMCParticle * particle);

//__________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  //-- get an instance of the algorithm factory

  LOG("Main",pINFO) << "Getting an instance of the AlgFactory";

  AlgFactory * algf = AlgFactory::Instance();

  //-- get the PythiaDecayer algorithm implementing the DecayModelI interface

  LOG("Main",pINFO)
          << "Asking the factory for the genie::PythiaDecayer\\Default alg.";

  const Algorithm * alg_base = algf->GetAlgorithm(
                                           "genie::PythiaDecayer","Default");

  const DecayModelI * pythia_decayer =
                                dynamic_cast<const DecayModelI *> (alg_base);

  //-- ask the algorithms for its name and its configuration parameter set

  LOG("Main",pINFO) << "Algorithm name = " << pythia_decayer->Id().Name();
  LOG("Main",pINFO) << "Parameter set  = " << pythia_decayer->Id().Config();


  //-- get the configuration registries and print them

  LOG("Main",pINFO) << "Getting/Printing the configuration registries";

  const Registry & conf_registry = pythia_decayer->GetConfig();

  LOG("Main",pINFO) << conf_registry;

  //-- Use the PythiaDecayer to decay a particle

  // Decaying particle 4-momentum

  TLorentzVector p4;

  p4.SetE(30);
  p4.SetTheta(0);
  p4.SetPhi(0);

  LOG("Main",pINFO) << "Decaying a pion with energy = " << p4.Energy();

  DecayerInputs_t dinp;

  dinp.PdgCode = kPdgPiP;
  dinp.P4      = &p4;

  TClonesArray * particle_list = pythia_decayer->Decay(dinp);

  // Print decay products

  LOG("Main",pINFO) << particle_list;  

  particle_list->Delete();

  delete particle_list;

  return 0;
}
//__________________________________________________________________________
ostream & operator<< (ostream & stream, const TClonesArray * particle_list)
{
  TMCParticle * p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (TMCParticle *) particle_iter.Next()) ) stream << p;

  return stream; 
}
//__________________________________________________________________________
ostream & operator<< (ostream & stream, const TMCParticle * particle)
{
  stream << endl
         << "name = " << particle->GetName()
         << " KF = "  << particle->GetKF() 
         << " KS = "  << particle->GetKS()
         << "(E = "   << particle->GetEnergy() 
         << ",Px = "  << particle->GetPx()
         << ",Py = "  << particle->GetPy()
         << ",Pz = "  << particle->GetPz() << ")";

  return stream;
}
//__________________________________________________________________________

