//____________________________________________________________________________
/*!

\program gtestDecay

\brief   test program used for testing / debugging the PYTHIA/JETSET
         interface and the GENIE DecayModelI algorithms

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
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
#include "Conventions/Units.h"
#include "Conventions/Constants.h"
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
  // Get the decayer
  LOG("Main",pINFO)
     << "Asking the AlgFactory for a genie::PythiaDecayer\\Default instance";
  AlgFactory * algf = AlgFactory::Instance();
  const DecayModelI * decayer =
             dynamic_cast<const DecayModelI *> (
                    algf->GetAlgorithm("genie::PythiaDecayer","Default"));

  // Decayer config print-out
  LOG("Main",pINFO) << "Algorithm name = " << decayer->Id().Name();
  LOG("Main",pINFO) << "Parameter set  = " << decayer->Id().Config();

  const Registry & conf_registry = decayer->GetConfig();
  LOG("Main",pINFO) << conf_registry;

  // Define a particle code/4-p and setup the decayer inputs

  int    pdgc  = kPdgPi0; // pi0
  double E     = 30.;     // GeV
  double theta = 0;  
  double phi   = 0;  

  DecayerInputs_t dinp;

  TLorentzVector p4;
  p4.SetE(E);
  p4.SetTheta(theta);
  p4.SetPhi(phi);

  dinp.PdgCode = pdgc;
  dinp.P4      = &p4;

  LOG("Main",pINFO) 
       << "Decaying a pion with energy = " << p4.Energy();

  // Perform the decay a few times & print-out the decay products
  const int ndec = 40;
  for(int idec = 0; idec < ndec; idec++) {

        LOG("Main",pINFO) << "*** Decay nu.: " << idec;

        // Decay
  	TClonesArray * particle_list = decayer->Decay(dinp);

        // Print decay products
        LOG("Main",pINFO) << particle_list;  

        // Clean-up
        particle_list->Delete();
        delete particle_list;
  }

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
         << "(E = "   << particle->GetEnergy() << " GeV"
         << ", Px = "  << particle->GetPx() << " GeV/c"
         << ", Py = "  << particle->GetPy() << " GeV/c"
         << ", Pz = "  << particle->GetPz() << " GeV/c) "
         << "(t = "   << particle->GetTime() /(units::mm) << " mm/c"
         << ", x = "  << particle->GetVx()   /(units::mm) << " mm"
         << ", y = "  << particle->GetVy()   /(units::mm) << " mm"
         << ", z = "  << particle->GetVz()   /(units::mm) << " mm)";

  return stream;
}
//__________________________________________________________________________

