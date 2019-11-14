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
#include <TIterator.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Physics/Hadronization/FragmRecUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
                       int pdg_code, const TClonesArray * const particle_list)
{
  int nparticles = 0;
  GHepParticle* p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (GHepParticle *) particle_iter.Next()) ) {
    if(p->Pdg() == pdg_code)  {
      if(p->Status()<10) nparticles++;
    }
  }
  return nparticles;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NParticles(
   int pdg_code, GHepStatus_t status, const TClonesArray * const particle_list)
{
  int nparticles = 0;
  GHepParticle* p = 0;

  TObjArrayIter particle_iter(particle_list);

  while( (p = (GHepParticle *) particle_iter.Next()) )
              if(p->Pdg() == pdg_code && p->Status() == status) nparticles++;

  return nparticles;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NPositives(const TClonesArray * const part_list)
{
// Find out the number of negative particles in the particle container

  TIter piter(part_list);

  GHepParticle * p = 0;
  int npos = 0;

  while( (p = (GHepParticle *) piter.Next()) )
         if( PDGLibrary::Instance()->Find(p->Pdg())->Charge() > 0 ) npos++;

  return npos;
}
//____________________________________________________________________________
int genie::utils::fragmrec::NNegatives(const TClonesArray * const part_list)
{
// Find out the number of negative particles in the particle container

  TIter piter(part_list);

  GHepParticle * p = 0;
  int nneg = 0;

  while( (p = (GHepParticle *) piter.Next()) )
         if( PDGLibrary::Instance()->Find(p->Pdg())->Charge() < 0 ) nneg++;

  return nneg;
}
//____________________________________________________________________________
void genie::utils::fragmrec::Print(const TClonesArray * const part_list)
{
  TIter piter(part_list);

  unsigned int i=0;
  GHepParticle * particle = 0;

  double sum_px = 0, sum_py = 0, sum_pz = 0, sum_E = 0;


  while( (particle = (GHepParticle *) piter.Next()) ) {

    sum_E  += (particle->Energy());
    sum_px += (particle->Px());
    sum_py += (particle->Py());
    sum_pz += (particle->Pz());

    SLOG("FragmRecUtils", pINFO)
        << "-> " << i++ << " " << particle->Name()
        << " PDG = " << particle->Pdg()
        << " status = " << particle->Status()
        << " moms = {"
        << particle->FirstMother() << ", " << particle->LastMother()
        << "} kids = {"
        << particle->FirstDaughter() << ", " << particle->LastDaughter()
        << "}(E = "  << particle->Energy()
        << ",Px = "  << particle->Px()
        << ",Py = "  << particle->Py()
        << ",Pz = "  << particle->Pz() << ")";
  }

  SLOG("FragmRecUtils", pINFO)
       << "SUMS: E = " << sum_E
       << ", px = " << sum_px << ", py = " << sum_py << ", pz = " << sum_pz;

}
//__________________________________________________________________________
