//____________________________________________________________________________
/*!

\namespace  genie::utils::fragmrec

\brief      Simple utilities for the Fragmentation Event Record.

            The Fragmentation event record is a TClonesArray of TMCParticles -
            equivalent to PYTHIA's PYJETS.

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    November 26, 2004

*/
//____________________________________________________________________________

#include <TMCParticle6.h>
#include <TIterator.h>

#include "Messenger/Messenger.h"
#include "Utils/FragmRecUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;

//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
                       int pdg_code, const TClonesArray * const particle_list)
{
  int nparticles = 0;

  TMCParticle * p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (TMCParticle *) particle_iter.Next()) )
                                     if(p->GetKF() == pdg_code)  nparticles++;

  return nparticles;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
           int pdg_code, int status, const TClonesArray * const particle_list)
{
  int nparticles = 0;

  TMCParticle * p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (TMCParticle *) particle_iter.Next()) )
              if(p->GetKF() == pdg_code && p->GetKS() == status) nparticles++;

  return nparticles;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NPositives(const TClonesArray * const part_list)
{
// Find out the number of negative particles in the particle container

  TIter piter(part_list);

  TMCParticle * p = 0;

  int npos = 0;

  while( (p = (TMCParticle *) piter.Next()) )
         if( PDGLibrary::Instance()->Find(p->GetKF())->Charge() > 0 ) npos++;

  return npos;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NNegatives(const TClonesArray * const part_list)
{
// Find out the number of negative particles in the particle container

  TIter piter(part_list);

  TMCParticle * p = 0;

  int nneg = 0;

  while( (p = (TMCParticle *) piter.Next()) )
         if( PDGLibrary::Instance()->Find(p->GetKF())->Charge() < 0 ) nneg++;

  return nneg;
}
//____________________________________________________________________________
void genie::utils::fragmrec::Print(const TClonesArray * const part_list)
{
  TIter piter(part_list);

  TMCParticle * particle = 0;

  double sum_px = 0, sum_py = 0, sum_pz = 0, sum_E = 0;

  while( (particle = (TMCParticle *) piter.Next()) ) {

    sum_E  += (particle->GetEnergy());
    sum_px += (particle->GetPx());
    sum_py += (particle->GetPy());
    sum_pz += (particle->GetPz());

    SLOG("FragmResUtils", pINFO)
            << " Name = " << particle->GetName()
            << " KF = "   << particle->GetKF()
            << " KS = "   << particle->GetKS()
            << "(E = "    << particle->GetEnergy()
            << ",Px = "   << particle->GetPx()
            << ",Py = "   << particle->GetPy()
            << ",Pz = "   << particle->GetPz() << ")";
  }

  SLOG("FragmResUtils", pINFO)
       << "SUMS: E = " << sum_E
       << ", px = " << sum_px << ", py = " << sum_py << ", pz = " << sum_pz;

}
//__________________________________________________________________________

