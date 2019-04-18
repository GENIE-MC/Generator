//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 26, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TIterator.h>

#include "Framework/Messenger/Messenger.h"
#include "Physics/Hadronization/FragmRecUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
                       int pdg_code, const TClonesArray * const particle_list)
{
  int nparticles = 0;
  TMCParticle* p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (TMCParticle *) particle_iter.Next()) ) {
    if(p->GetKF() == pdg_code)  {
      if(p->GetKS()<10) nparticles++;
    }
  }
  return nparticles;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
           int pdg_code, int status, const TClonesArray * const particle_list)
{
  int nparticles = 0;
  TMCParticle* p = 0;

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

  unsigned int i=0;
  TMCParticle * particle = 0;

  double sum_px = 0, sum_py = 0, sum_pz = 0, sum_E = 0;

  
  while( (particle = (TMCParticle *) piter.Next()) ) {

    sum_E  += (particle->GetEnergy());
    sum_px += (particle->GetPx());
    sum_py += (particle->GetPy());
    sum_pz += (particle->GetPz());

    SLOG("FragmRecUtils", pINFO)
        << "-> " << i++ << " " << particle->GetName()
        << " KF = " << particle->GetKF()
        << " KS = " << particle->GetKS()
        << " mom = " << particle->GetParent()
        << " kids = {" 
        << particle->GetFirstChild() << ", " << particle->GetLastChild() 
        << "}(E = "  << particle->GetEnergy()
        << ",Px = "  << particle->GetPx()
        << ",Py = "  << particle->GetPy()
        << ",Pz = "  << particle->GetPz() << ")";
  }

  SLOG("FragmRecUtils", pINFO)
       << "SUMS: E = " << sum_E
       << ", px = " << sum_px << ", py = " << sum_py << ", pz = " << sum_pz;

}
//__________________________________________________________________________

