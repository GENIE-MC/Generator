//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool - November 19, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TLorentzVector.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/HadronTransport/NucBindEnergyAggregator.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/MathUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucBindEnergyAggregator::NucBindEnergyAggregator() :
EventRecordVisitorI("genie::NucBindEnergyAggregator")
{

}
//___________________________________________________________________________
NucBindEnergyAggregator::NucBindEnergyAggregator(string config) :
EventRecordVisitorI("genie::NucBindEnergyAggregator", config)
{

}
//___________________________________________________________________________
NucBindEnergyAggregator::~NucBindEnergyAggregator()
{

}
//___________________________________________________________________________
void NucBindEnergyAggregator::ProcessEventRecord(GHepRecord * evrec) const
{
  // Return if the neutrino was not scatterred off a nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("Nuclear", pINFO)
        << "No nuclear target found - Not subtracting any binding energy";
    return;
  }

  // Loop over particles, find final state nucleons for which a removal 
  // energy has been set and subtract it from their energy
  TIter stdhep_iter(evrec);
  GHepParticle * p = 0;

  while( (p = (GHepParticle * ) stdhep_iter.Next()) ) {

     bool is_nucleon   = pdg::IsNeutronOrProton(p->Pdg());  
     bool in_fin_state = (p->Status() == kIStStableFinalState);
     bool had_bind_e   = (p->RemovalEnergy() > 0.);

     bool handle = is_nucleon && in_fin_state && had_bind_e;
     if(!handle) continue;

     //-- ask for the binding energy set by the nuclear model
     double bindE = p->RemovalEnergy();
     LOG("Nuclear", pINFO) << "Binding energy = " << bindE;

     //-- subtract this energy from the final state nucleon
     LOG("Nuclear", pINFO)
          << "Subtracting the binding energy from the escaped nucleon";

     double M  = p->Mass();
     double En = p->Energy();
     double KE = En-M;

     LOG("Nuclear", pINFO)  << "Kinetic energy before subtraction = " << KE;
     KE -= bindE;

     LOG("Nuclear", pINFO) << "Kinetic energy after subtraction = " << KE;

     En = KE+M;

     if(En>M || !fAllowRecombination) { 
         double pmag_old = p->P4()->P();
         double pmag_new = TMath::Sqrt(utils::math::NonNegative(En*En-M*M));
         double scale    = pmag_new / pmag_old;
         LOG("Nuclear", pINFO) 
             << "|pnew| = " << pmag_new << ", |pold| = " << pmag_old
             << ", scale = " << scale;

         double pxn = scale * p->Px();
         double pyn = scale * p->Py();
         double pzn = scale * p->Pz();

         double pxb = (1-scale) * p->Px();
         double pyb = (1-scale) * p->Py();
         double pzb = (1-scale) * p->Pz();

          p->SetEnergy ( En  );
          p->SetPx     ( pxn );
          p->SetPy     ( pyn );
          p->SetPz     ( pzn );
          p->SetRemovalEnergy(0);

          //-- and add a GHEP entry to record this in the event record and
          //   conserve energy/momentum
          LOG("Nuclear", pINFO)
             << "Adding a [BindingE] to account for nuclear binding energy";

          evrec->AddParticle(kPdgBindino, kIStStableFinalState, 
      		                   -1,-1,-1,-1, pxb,pyb,pzb,bindE, 0,0,0,0);
     } else {
          LOG("Nuclear", pNOTICE)
               << "Nucleon is above the Fermi sea but can't escape the nucleus";
          LOG("Nuclear", pNOTICE)
               << "Recombining remnant nucleus + f/s nucleon";


          // find the remnant nucleus if corresponding cascade modules provide it
          int rnucpos = evrec->RemnantNucleusPosition();

          if (rnucpos != -1)
          {
            GHepParticle * rnucl = evrec->Particle(rnucpos);

            // mark both the remnant nucleus and the final state nucleon as 
            // intermediate states
            rnucl -> SetStatus(kIStIntermediateState);
            p     -> SetStatus(kIStIntermediateState);

            // figure out the recombined nucleus PDG code
            int Z = rnucl->Z();
            int A = rnucl->A();
            if(pdg::IsProton(p->Pdg())) Z++;
            A++;
            int ipdgc = pdg::IonPdgCode(A,Z);

            // add-up their 4-momenta
            double pxnuc = rnucl->Px() + p->Px();
            double pynuc = rnucl->Py() + p->Py();
            double pznuc = rnucl->Pz() + p->Pz();
            double Enuc  = rnucl->E()  + p->E();

            evrec->AddParticle(ipdgc, kIStStableFinalState,
                        rnucpos,-1,-1,-1, pxnuc,pynuc,pznuc,Enuc, 0,0,0,0);
          }
          else
          {
              // Find hadronic blob which is provided by Intranuke models
              GHepParticle *hblob = evrec->FindParticle(kPdgHadronicBlob  , kIStFinalStateNuclearRemnant, 0);
              // If there is no hadronic blob do nothing
              if (hblob == nullptr) return;
              p     -> SetStatus(kIStIntermediateState);
              
              // add-up their 4-momenta
              double pxhblob = hblob->Px() + p->Px();
              double pyhblob = hblob->Py() + p->Py();
              double pzhblob = hblob->Pz() + p->Pz();
              double Ehblob  = hblob->E()  + p->E();
              hblob->SetMomentum(pxhblob, pyhblob, pzhblob, Ehblob);
          }
     }
  }
}
//___________________________________________________________________________
/*
GHepParticle * NucBindEnergyAggregator::FindMotherNucleus(
                                         int ipos, GHepRecord * evrec) const
{
  GHepParticle * p = evrec->Particle(ipos);

  //-- get its mothet
  int mother_pos = p->FirstMother();

  //-- if mother is set
  if(mother_pos != -1) {
     GHepParticle * mother = evrec->Particle(mother_pos);

     //-- check its status
     if( mother->Status() == kIStNucleonTarget ) {

        //-- get the mother's mother
        int grandmother_pos = mother->FirstMother();

        //-- if grandmother is set get its PDG code a check if it is an ion
        if(grandmother_pos != -1) {
             GHepParticle * grandmother =
                                     evrec->Particle(grandmother_pos);

             int grandmother_pdgc = grandmother->Pdg();
             if( pdg::IsIon(grandmother_pdgc) ) return grandmother;

        } // gmpos != -1
     } // mother-status
  }  //mpos != -1

  return 0;
} */
//___________________________________________________________________________
void NucBindEnergyAggregator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucBindEnergyAggregator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucBindEnergyAggregator::LoadConfig(void)
{
	GetParamDef("AllowNuclRecombination",fAllowRecombination, true  ) ;

}
//____________________________________________________________________________

