//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models.

          Is a concrete implementation of the HadronizationModelI interface.
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

*/
//____________________________________________________________________________

#include <TMCParticle6.h>

#include "Fragmentation/PythiaHadronization.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//-- the actual PYTHIA call

extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
PythiaHadronization::PythiaHadronization() :
HadronizationModelI("genie::PythiaHadronization")
{
  fPythia = new TPythia6();  
}
//____________________________________________________________________________
PythiaHadronization::PythiaHadronization(string config) :
HadronizationModelI("genie::PythiaHadronization", config)
{
  fPythia = new TPythia6();  
}
//____________________________________________________________________________
PythiaHadronization::~PythiaHadronization() 
{ 

}
//____________________________________________________________________________
void PythiaHadronization::Initialize(void) const
{

}  
//____________________________________________________________________________
TClonesArray * PythiaHadronization::Hadronize(
                                        const Interaction * interaction) const
{
  const Kinematics & kinematics = interaction->GetKinematics();
  double W = kinematics.W();

  const InitialState & init_state = interaction->GetInitialState();

  int hit_nucleon = init_state.GetTarget().StruckNucleonPDGCode();
  int hit_quark   = init_state.GetTarget().StruckQuarkPDGCode();
  int diquark     = 0;

  //-- check for hit-quark assignment
  
  if(hit_quark==0) {

     // no hit-quark assignement - selecting a random one for the input nucleon
     // (should actually use PDFs here...)
     RandomGen * rnd = RandomGen::Instance();

     double x = rnd->Random1().Rndm();
     if (pdg::IsProton(hit_nucleon)) {
         hit_quark = kPdgUQuark;
         if(x < 0.3333) hit_quark = kPdgDQuark;
     } else if(pdg::IsNeutron(hit_nucleon)) {
         hit_quark = kPdgUQuark;
         if(x < 0.6666) hit_quark = kPdgDQuark;
     }
  }

  //-- PYTHIA->HADRONIZE:

  int ip = 0;

  //-- hit nucleon = proton

  if ( pdg::IsProton(hit_nucleon) ) {

       /* u + ud */
       if ( pdg::IsUQuark(hit_quark) )
       {
          diquark = kPdgUDDiquarkS1;
          py2ent_(&ip, &hit_quark, &diquark, &W);
       }
       /* d + uu */
       else if ( pdg::IsDQuark(hit_quark) )
       {
          diquark = kPdgUUDiquarkS1;
          py2ent_(&ip, &hit_quark, &diquark, &W);
       }
       else {
         LOG("PythiaHad", pERROR) << "Can not handle quark: " << hit_quark;
       }
  }

  //-- hit nucleon = neutron

  else if ( pdg::IsNeutron(hit_nucleon) ) {

       /* u + dd */
       if ( pdg::IsUQuark(hit_quark) )
       {
          diquark = kPdgDDDiquarkS1;
          py2ent_(&ip, &hit_quark, &diquark, &W);
       }          
        /* d + ud */
       else if ( pdg::IsDQuark(hit_quark) )
       {
          diquark = kPdgUDDiquarkS1;
          py2ent_(&ip, &hit_quark, &diquark, &W);
       }       
       else {
          LOG("PythiaHad", pERROR) << "Can not handle quark: " << hit_quark;
       }

  } else  {    
      LOG("PythiaHad", pERROR) << "Can not handle nucleon: " << hit_nucleon;
  }

  fPythia->GetPrimaries();

  TClonesArray * pythia_particles =
                           (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  TClonesArray * particle_list = new TClonesArray(
                             "TMCParticle", pythia_particles->GetEntries() );
                             
  register unsigned int i = 0;
  TMCParticle * particle = 0;
  
  TIter particle_iter(pythia_particles);

  while( (particle = (TMCParticle *) particle_iter.Next()) ) {
    
       LOG("PythiaHad", pINFO)
               << "Adding final state particle PDGC = " << particle->GetKF();

       new ( (*particle_list)[i++] ) TMCParticle(*particle);
  }

  particle_list->SetOwner(true);
    
  return particle_list;
}
//____________________________________________________________________________
