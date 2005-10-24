//____________________________________________________________________________
/*!

\class   genie::UnstableParticleDecayer

\brief   Decays unstable particles found in the generated event record.

         After the interaction vertex generation it visits the event record
         and it decays the unstable particles using an externally specified
         particle decay model. The decay products are added to the event
         record and the status of parent particle is toggled. \n

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 17, 2004

*/
//____________________________________________________________________________

#include <algorithm>

#include <TParticlePDG.h>
#include <TMCParticle6.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Decay/DecayModelI.h"
#include "EVGModules/UnstableParticleDecayer.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using std::count;

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer() :
EventRecordVisitorI("genie::UnstableParticleDecayer")
{

}
//___________________________________________________________________________
UnstableParticleDecayer::UnstableParticleDecayer(string config) :
EventRecordVisitorI("genie::UnstableParticleDecayer", config)
{

}
//___________________________________________________________________________
UnstableParticleDecayer::~UnstableParticleDecayer()
{

}
//___________________________________________________________________________
void UnstableParticleDecayer::ProcessEventRecord(GHepRecord * evrec) const
{
  //-- Get the specified decay model

  const DecayModelI * decayer = dynamic_cast<const DecayModelI *>
                      (this->SubAlg("decayer-alg-name","decayer-param-set"));

  //-- Loop over particles, find unstable ones and decay them

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  unsigned int ipos = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {

     if( this->ToBeDecayed(p) ) {
        LOG("ParticleDecayer", pINFO)
              << "\n Attempt to decay unstable particle: " << p->GetName();

        //-- Get the parent particle 4-momentum

        TLorentzVector p4(p->Px(), p->Py(), p->Pz(), p->E());

        //-- Decay it & retrieve the decay products
        //   The decayer might not be able to handle it - in which case it
        //   should return a NULL decay products container

        DecayerInputs_t dinp;

        dinp.PdgCode = p->PdgCode();
        dinp.P4      = &p4;

        TClonesArray * decay_products = decayer->Decay(dinp);

        //-- Check whether the particle was decayed
        if(decay_products) {
           LOG("ParticleDecayer", pINFO) << "The particle was decayed";

           //-- Mark it as a 'decayed state' & add its daughter links
           p->SetStatus(kIstDecayedState);

           //-- Loop over the daughter and add them to the event record
           this->CopyToEventRecord(decay_products, evrec, ipos);

           decay_products->Delete();
           delete decay_products;
        }// !=0
     }// to be decayed?
     ipos++;

  } // loop over particles
}
//___________________________________________________________________________
bool UnstableParticleDecayer::ToBeDecayed(GHepParticle * particle) const
{
   if( particle->PdgCode() != 0 &&
                particle->Status() == kIStStableFinalState)
                                           return this->IsUnstable(particle);
   return false;
}
//___________________________________________________________________________
bool UnstableParticleDecayer::IsUnstable(GHepParticle * particle) const
{
  int pdg_code = particle->PdgCode();

  TParticlePDG * ppdg = PDGLibrary::Instance()->Find(pdg_code);

  //-- Get the specified maximum lifetime tmax (decay with lifetime < tmax)

  double tmax /* sec */ =
          (fConfig->Exists("max-lifetime-for-unstables")) ?
                    fConfig->GetDouble("max-lifetime-for-unstables") : 1e-10;

   if( ppdg->Lifetime() < tmax ) { /*return true*/ };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - temporary/
  // ROOT's TParticlepdg::Lifetim e() does not work properly
  // do something else instead (temporarily)

  int particles_to_decay[] = {
          kPdgPi0, kPdgLambdacP, kPdgSigmacP, kPdgSigmacPP, kPdgDPlus,
          kPdgDMinus, kPdgD0, kPdgD0Bar, kPdgDsPlus, kPdgDsMinus};

  const int N = sizeof(particles_to_decay) / sizeof(int);

  int matches = count(particles_to_decay, particles_to_decay+N, pdg_code);

  if(matches > 0 || utils::res::IsBaryonResonance(pdg_code)) return true;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - /temporary

  return false;
}
//___________________________________________________________________________
void UnstableParticleDecayer::CopyToEventRecord(TClonesArray *
                decay_products, GHepRecord * evrec, int mother_pos) const
{
  TMCParticle * dpmc = 0;

  TObjArrayIter decay_iter(decay_products);

  while( (dpmc = (TMCParticle *) decay_iter.Next()) ) {

     TLorentzVector vdummy(0,0,0,0); // position 4-vector

     double px = dpmc->GetPx();
     double py = dpmc->GetPy();
     double pz = dpmc->GetPz();
     double E  = dpmc->GetEnergy();

     TLorentzVector p4(px, py, pz, E); // momentum 4-vector

     int          pdg    = dpmc->GetKF();
     GHepStatus_t status = GHepStatus_t (dpmc->GetKS());

     //-- Only add the decay products - the mother particle already exists
     //   (mother's daughter list will be automatically updated with every
     //    AddParticle() call)
     if(status == kIStStableFinalState) {
         evrec->AddParticle(pdg, status, mother_pos,-1,-1,-1, p4, vdummy);
     }
  }
}
//___________________________________________________________________________
