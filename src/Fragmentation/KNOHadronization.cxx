//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TMCParticle6.h>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Decay/DecayModelI.h"
#include "Fragmentation/KNOHadronization.h"
#include "Fragmentation/MultiplicityProbModelI.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils::print;

//____________________________________________________________________________
KNOHadronization::KNOHadronization() :
HadronizationModelI("genie::KNOHadronization")
{

}
//____________________________________________________________________________
KNOHadronization::KNOHadronization(string config) :
HadronizationModelI("genie::KNOHadronization", config)
{

}
//____________________________________________________________________________
KNOHadronization::~KNOHadronization()
{

}
//____________________________________________________________________________
void KNOHadronization::Initialize(void) const
{

}
//____________________________________________________________________________
TClonesArray * KNOHadronization::Hadronize(
                                        const Interaction * interaction) const
{
// Generate the hadronic system in a neutrino interaction using a KNO-based 
// model. 

  RandomGen * rnd = RandomGen::Instance();
  vector<int> * pdgcv   = 0;
  double * mass         = 0;
  unsigned int min_mult = 2;
  unsigned int mult     = 0;

  //----- Available invariant mass
  double W = interaction->GetKinematics().W();
  LOG("KNOHad", pINFO) << "W = " << W << " GeV";

  if(W <= kNucleonMass+kPionMass) {
     LOG("KNOHad", pWARN) 
        << "Low invariant mass, W = " << W << " GeV! Returning a null list";
     return 0;
  }

  //----- Init event weight (to be set if producing weighted events)
  fWeight = 1.;

  //----- Get the charge that the hadron shower needs to have so as to
  //      conserve charge in the interaction

  int maxQ = this->HadronShowerCharge(interaction);
  LOG("KNOHad", pINFO) << "Hadron Shower Charge = " << maxQ;

 //----- Build the multiplicity probabilities for the input interaction

  LOG("KNOHad", pDEBUG) << "Building Multiplicity Probability distribution";
  LOG("KNOHad", pDEBUG) << *interaction;
  const TH1D & mprob = fMultProbModel->ProbabilityDistribution(interaction);

  if(mprob.Integral("width")<=0) {
    LOG("KNOHad", pERROR) << "Empty multiplicity probability distribution!";
    return 0;
  }

  //----- FIND AN ALLOWED SOLUTION FOR THE HADRONIC FINAL STATE

  TLorentzVector p4(0,0,0,W);

  bool allowed_state=false;
  register unsigned int itry = 0;

  while(!allowed_state) 
  {
    itry++;

    //-- Go in error if a solution has not been found after many attempts
    if(itry>kMaxKNOHadSystIterations) {
       LOG("KNOHad", pERROR) 
         << "Couldn't generate hadronic multiplicities after: " 
                                               << itry << " attempts!";
       return 0;
    }

    //-- Generate a hadronic multiplicity 
    mult = TMath::Nint( mprob.GetRandom() );

    LOG("KNOHad", pINFO) << "Hadron multiplicity  = " << mult;

    //-- Check that the generated multiplicity is consistent with the charge
    //   thet the hadronic shower is required to have - else retry
    if(mult < (unsigned int) TMath::Abs(maxQ)) {
       LOG("KNOHad", pWARN) 
        << "Multiplicity not enough to generate hadronic charge! Retrying.";
      allowed_state = false;
      continue;
    }

    //-- Force a min multiplicity
    //   This should never happen if the multiplicity probability distribution
    //   was properly built
    if(mult < min_mult) {
      if(fForceMinMult) {
        LOG("KNOHad", pWARN) 
           << "Low generated multiplicity: " << mult 
              << ". Forcing to minimum accepted multiplicity: " << min_mult;
       mult = min_mult;
      } else {
        LOG("KNOHad", pFATAL) 
           << "Generated multiplicity: " << mult << " is too low - Quitting";
        return 0;
      }
    }

    //-- Determine what kind of particles we have in the final state
    pdgcv = this->GenerateFSHadronCodes(mult, maxQ, W);

    LOG("KNOHad", pNOTICE) 
         << "Generated multiplicity (@ W = " << W << "): " << pdgcv->size();

    // muliplicity might have been forced to smaller value if the invariant
    // mass of the hadronic system was not sufficient
    mult = pdgcv->size(); // update for potential change

    //-- Set requested decay to a phase space generator
    vector<int>::const_iterator pdg_iter;
    mass = new double[pdgcv->size()];
    int i = 0;
    double msum=0;
    LOG("KNOHad", pDEBUG) << "Current hadronic system:";
    for(pdg_iter = pdgcv->begin(); pdg_iter != pdgcv->end(); ++pdg_iter) {
      int pdgc = *pdg_iter;
      double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
      msum += m;
      mass[i++] = m;
      LOG("KNOHad", pDEBUG) << "- PDGC=" << pdgc << ", m=" << m << " GeV";
    }
    bool permitted = fPhaseSpaceGenerator.SetDecay(p4, mult, mass);

    //-- Check that the requested decay is permitted
    if(!permitted) {
       LOG("KNOHad", pWARN) << "*** Decay forbidden by kinematics! ***";
       LOG("KNOHad", pWARN) << "sum{mass} = " << msum << ", W = " << W;
       LOG("KNOHad", pWARN) << "Discarding hadronic system & re-trying!";
       delete pdgcv;
       delete [] mass;
       allowed_state = false;
       continue;
    }

    allowed_state = true;
    LOG("KNOHad", pNOTICE) 
             << "Found an allowed hadronic state @ W=" << W 
                                              << " multiplicity=" << mult;
  } // attempts

  //----- DECAY HADRONIC FINAL STATE

  //-- get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int i=0; i<200; i++) {
     double w = fPhaseSpaceGenerator.Generate();   
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);

  LOG("KNOHad", pNOTICE) 
     << "Max phase space gen. weight @ current hadronic system: " << wmax;

  if(fGenerateWeighted) 
  {
    // *** generating weighted decays ***
    double w = fPhaseSpaceGenerator.Generate();   
    fWeight *= TMath::Max(w/wmax, 1.);
  }
  else 
  {
    // *** generating un-weighted decays ***
     wmax *= 1.2;
     bool accept_decay=false;
     itry=0;

     while(!accept_decay) 
     {
       itry++;
       assert(itry<kMaxUnweightDecayIterations);

       double w  = fPhaseSpaceGenerator.Generate();   
       double gw = wmax * rnd->RndHadro().Rndm();

       LOG("KNOHad", pINFO) << "Decay weight = " << w << " / R = " << gw;

       accept_decay = (gw<=w);
     }
  }

  //-- insert final state products into a TClonesArray of TMCParticles
  //   and return it
  LOG("KNOHad", pDEBUG)
       << "Creating output particle list TClonesArray of TMCParticles";

  TClonesArray * particle_list = new TClonesArray("TMCParticle", mult);

  for(unsigned int i = 0; i < pdgcv->size(); i++) {

     int pdgc = (*pdgcv)[i];

     //-- get the 4-momentum of the i-th final state particle
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(i);

     //-- push TMCParticle in the particle list (TClonesArray)
     LOG("KNOHad", pDEBUG)
            << "Adding final state particle PDGC = " << pdgc
               << " with mass = " << mass[i] << " GeV";

     new ( (*particle_list)[i] ) TMCParticle
        (
           1,               /* KS Code                          */
           pdgc,            /* PDG Code                         */
           0,               /* parent particle                  */
           0,               /* first child particle             */
           0,               /* last child particle              */
           p4fin->Px(),     /* 4-momentum: px component         */
           p4fin->Py(),     /* 4-momentum: py component         */
           p4fin->Pz(),     /* 4-momentum: pz component         */
           p4fin->Energy(), /* 4-momentum: E  component         */
           mass[i],         /* particle mass                    */
           0,               /* production vertex 4-vector: vx   */
           0,               /* production vertex 4-vector: vy   */
           0,               /* production vertex 4-vector: vz   */
           0,               /* production vertex 4-vector: time */
           0                /* lifetime                         */
        );
  }

  // handle unstable particle decays (if requested)
  this->HandleDecays(particle_list);

  // the container 'owns' its elements
  particle_list->SetOwner(true);

  delete pdgcv;
  delete [] mass;

  return particle_list;
}
//____________________________________________________________________________
double KNOHadronization::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
void KNOHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KNOHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KNOHadronization::LoadConfig(void)
{
// Read configuration options or set defaults

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fMultProbModel = 0;
  fDecayer       = 0;

  // Force decays of unstable hadronization products?
  fForceDecays  = fConfig->GetBoolDef("force-decays", false);

  // Force minimum multiplicity (if generated less than that) or abort?
  fForceMinMult = fConfig->GetBoolDef("force-min-multiplicity", true);

  // Probabilities for producing hadron pairs

  //-- pi0 pi0
  fPpi0 = fConfig->GetDoubleDef(
                         "prob-fs-pi0-pair", gc->GetDouble("KNO-ProbPi0Pi0")); 
  //-- pi+ pi-
  fPpic = fConfig->GetDoubleDef(
            "prob-fs-piplus-piminus", gc->GetDouble("KNO-ProbPiplusPiminus")); 

  //-- K+  K-
  fPKc  = fConfig->GetDoubleDef(
                "prob-fs-Kplus-Kminus", gc->GetDouble("KNO-ProbKplusKminus")); 

  //-- K0 K0bar
  fPK0  = fConfig->GetDoubleDef(
                        "prob-fs-K0-K0bar", gc->GetDouble("KNO-ProbK0K0bar")); 

  // Multiplicity probability model
  fMultProbModel = dynamic_cast<const MultiplicityProbModelI *> (
    this->SubAlg("multiplicity-prob-alg-name", "multiplicity-prob-param-set"));
  assert(fMultProbModel);

  // Decay unstable particles now or leave it for later? Which decayer to use?
  if(fForceDecays) {
      fDecayer = dynamic_cast<const DecayModelI *> (
                       this->SubAlg("decayer-alg-name", "decayer-param-set"));
      assert(fDecayer);
  }

  // Generated weighted or un-weighted hadronic systems
  fGenerateWeighted = fConfig->GetBoolDef("generate-weighted", false);
}
//____________________________________________________________________________
vector<int> * KNOHadronization::GenerateFSHadronCodes(
                                   int multiplicity, int maxQ, double W) const
{
  // Code translated & adapted from NeuGEN. Original author: H.Gallagher
  // Non-trival changes -- needs validation against original code & data
  // Generate final state hadrons

  // PDG Library
  PDGLibrary * pdg = PDGLibrary::Instance();

  // vector to add final state hadron PDG codes
  vector<int> * pdgc = new vector<int>;
  pdgc->reserve(multiplicity);
  int hadrons_to_add = multiplicity;

  //----- Assign baryon as p or n.
  pdgc->push_back( GenerateBaryonPdgCode(multiplicity, maxQ) );

  // update number of hadrons to add, available shower charge & invariant mass
  if( (*pdgc)[0] == kPdgProton ) maxQ -= 1;
  hadrons_to_add--;
  W -= pdg->Find( (*pdgc)[0] )->Mass();

  //----- Assign remaining hadrons up to n = multiplicity

  //-- Handle charge imbalance
  while(maxQ != 0) {

     if (maxQ < 0) {
        // Need more negative charge
        LOG("KNOHad", pDEBUG) << "Need more negative charge -> Adding a pi-";
        pdgc->push_back( kPdgPiMinus );

        // update n-of-hadrons to add, avail. shower charge & invariant mass
        maxQ += 1;
        hadrons_to_add--;

        W -= pdg->Find(kPdgPiMinus)->Mass();

     } else if (maxQ > 0) {
        // Need more positive charge
        LOG("KNOHad", pDEBUG) << "Need more positive charge -> Adding a pi+";
        pdgc->push_back( kPdgPiPlus );

        // update n-of-hadrons to add, avail. shower charge & invariant mass
        maxQ -= 1;
        hadrons_to_add--;

        W -= pdg->Find(kPdgPiPlus)->Mass();
     }
  }

  //-- Add remaining neutrals or pairs up to the generated multiplicity
  if(maxQ == 0) {

     LOG("KNOHad", pDEBUG) 
       << "Hadronic charge balanced. Now adding only neutrals or +- pairs";

     // Final state has correct charge.
     // Now add pi0 or pairs (pi0 pi0 / pi+ pi- / K+ K- / K0 K0bar) only

     // Masses of particle pairs
     double M2pi0 = 2 * pdg -> Find (kPdgPi0    ) -> Mass();
     double M2pic =     pdg -> Find (kPdgPiPlus ) -> Mass() +
                        pdg -> Find (kPdgPiMinus) -> Mass();
     double M2Kc  =     pdg -> Find (kPdgKPlus  ) -> Mass() +
                        pdg -> Find (kPdgKMinus ) -> Mass();
     double M2K0  = 2 * pdg -> Find (kPdgK0     ) -> Mass();

     // Prevent multiplicity overflow.
     // Check if we have an odd number of hadrons to add.
     // If yes, add a single pi0 and then go on and add pairs

     if( hadrons_to_add > 0 && hadrons_to_add % 2 == 1 ) {

        LOG("KNOHad", pDEBUG)
                  << "Odd number of hadrons left to add -> Adding a pi0";
        pdgc->push_back( kPdgPi0 );

        // update n-of-hadrons to add & available invariant mass
        hadrons_to_add--;
        W -= pdg->Find(kPdgPi0)->Mass();
     }

     // Now add pairs (pi0 pi0 / pi+ pi- / K+ K- / K0 K0bar)

     assert( hadrons_to_add % 2 == 0 ); // even number

     RandomGen * rnd = RandomGen::Instance();

     while(hadrons_to_add > 0 && W >= M2pi0) {

         double x = rnd->RndHadro().Rndm();
         LOG("KNOHad", pDEBUG) << "rndm = " << x;

         if (x >= 0 && x < fPpi0) {
            //-------------------------------------------------------------
            // Add a pi0-pair
            LOG("KNOHad", pDEBUG) << " -> Adding a pi0pi0 pair";

            pdgc->push_back( kPdgPi0 );
            pdgc->push_back( kPdgPi0 );

            hadrons_to_add -= 2; // update the number of hadrons to add
            W -= M2pi0; // update the available invariant mass
            //-------------------------------------------------------------

         } else if (x < fPpi0 + fPpic) {
            //-------------------------------------------------------------
            // Add a piplus - piminus pair if there is enough W
            if(W >= M2pic) {

                LOG("KNOHad", pDEBUG) << " -> Adding a pi+pi- pair";

                pdgc->push_back( kPdgPiPlus  );
                pdgc->push_back( kPdgPiMinus );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2pic; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) 
                  << "Not enough mass for a pi+pi-: trying something else";
            }
            //-------------------------------------------------------------

         } else if (x < fPpi0 + fPpic + fPKc) {
            //-------------------------------------------------------------
            // Add a Kplus - Kminus pair if there is enough W
            if(W >= M2Kc) {

                LOG("KNOHad", pDEBUG) << " -> Adding a K+K- pair";

                pdgc->push_back( kPdgKPlus  );
                pdgc->push_back( kPdgKMinus );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2Kc; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) 
                     << "Not enough mass for a K+K-: trying something else";
            }
            //-------------------------------------------------------------

         } else if (x <= fPpi0 + fPpic + fPKc + fPK0) {
            //-------------------------------------------------------------
            // Add a K0 - K0bar pair if there is enough W
            if( W >= M2K0 ) {
                LOG("KNOHad", pDEBUG) << " -> Adding a K0 K0bar pair";

                pdgc->push_back( kPdgK0 );
                pdgc->push_back( kPdgK0 );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2K0; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) 
                 << "Not enough mass for a K0 K0bar: trying something else";
            }
            //-------------------------------------------------------------

         } else {
            LOG("KNOHad", pERROR)
                 << "Hadron Assignment Probabilities do not add up to 1!!";
            exit(1);
         }

         // make sure it has enough invariant mass to reach the
         // given multiplicity, even by adding only the lightest
         // hadron pairs (pi0's)
         // Otherwise force a lower multiplicity.
         if(W < M2pi0) hadrons_to_add = 0;

     } // while there are more hadrons to add
  } // if charge is balanced (maxQ == 0)

  return pdgc;
}
//____________________________________________________________________________
int KNOHadronization::GenerateBaryonPdgCode(int multiplicity, int maxQ) const
{
  // Code translated & adapted from NeuGEN. Original author: H.Gallagher
  // Assign baryon as p or n.
  // Force it for ++ and - I=3/2 at mult. = 2

  RandomGen * rnd = RandomGen::Instance();

  // initialize to neutron & then change it to proton if you must
  int pdgc = kPdgNeutron;
  double x = rnd->RndHadro().Rndm();

  //-- available hadronic system charge = 2
  if(maxQ == 2) {
     // for multiplicity == 2, force it to p
     if(multiplicity == 2) pdgc = kPdgProton;
     else {
           if(x < 0.66667) pdgc = kPdgProton;
     }
  }

  //-- available hadronic system charge = 1
  if(maxQ == 1) {
     if(multiplicity == 2) {
              if(x < 0.33333) pdgc = kPdgProton;
     } else {
              if(x < 0.50000) pdgc = kPdgProton;
     }
  }

  //-- available hadronic system charge = 0
  if(maxQ == 0) {
     if(multiplicity == 2) {
              if(x < 0.66667) pdgc = kPdgProton;
     } else {
              if(x < 0.50000) pdgc = kPdgProton;
     }
  }

  //-- available hadronic system charge = -1
  if(maxQ == -1) {
     // for multiplicity == 2, force it to n
     if(multiplicity != 2) {
           if(x < 0.33333) pdgc = kPdgProton;
     }
  }

  LOG("KNOHad", pDEBUG) 
       << " -> Adding a " << (( pdgc == kPdgProton ) ? "proton" : "neutron");

  return pdgc;
}
//____________________________________________________________________________
int KNOHadronization::HadronShowerCharge(const Interaction* interaction) const
{
// Returns the hadron shower charge in units of +e
// HadronShowerCharge = Q{initial} - Q{final state primary lepton}
// eg in v p -> l- X the hadron shower charge is +2
//    in v n -> l- X the hadron shower charge is +1
//    in v n -> v  X the hadron shower charge is  0
//
  int HadronShowerCharge = 0;

  // find out the charge of the final state lepton
  double ql = interaction->GetFSPrimaryLepton()->Charge() / 3.;

  // get the initial state, ask for the hit-nucleon and get
  // its charge ( = initial state charge for vN interactions)
  const InitialState & init_state = interaction->GetInitialState();
  int hit_nucleon = init_state.GetTarget().StruckNucleonPDGCode();

  assert( pdg::IsProton(hit_nucleon) || pdg::IsNeutron(hit_nucleon) );

  // Ask PDGLibrary for the nucleon charge
  double qinit = PDGLibrary::Instance()->Find(hit_nucleon)->Charge() / 3.;

  // calculate the hadron shower charge
  HadronShowerCharge = (int) ( qinit - ql );

  return HadronShowerCharge;
}
//____________________________________________________________________________
void KNOHadronization::HandleDecays(TClonesArray * plist) const
{
// Handle decays of unstable particles if requested through the XML config.
// The default is not to decay the particles at this stage (during event
// generation, the UnstableParticleDecayer event record visitor decays what
// is needed to be decayed later on). But, when comparing various models
// (eg PYTHIA vs KNO) independently and not within the full MC simulation
// framework it might be necessary to force the decays at this point.

  if (fForceDecays) {
     assert(fDecayer);

     //-- loop through the fragmentation event record & decay unstables
     int idecaying   = -1; // position of decaying particle
     TMCParticle * p =  0; // current particle

     TIter piter(plist);
     while ( (p = (TMCParticle *) piter.Next()) ) {

        idecaying++;
        int status = p->GetKS();

        // bother for final state particle only
        if(status < 10) {

          // until ROOT's T(MC)Particle(PDG) Lifetime() is fixed, decay only
          // pi^0's
          if ( p->GetKF() == kPdgPi0 ) {

               DecayerInputs_t dinp;

               TLorentzVector p4;
               p4.SetPxPyPzE(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());

               dinp.PdgCode = p->GetKF();
               dinp.P4      = &p4;

               TClonesArray * decay_products = fDecayer->Decay(dinp);

               if(decay_products) {
                  //--  mark the parent particle as decayed & set daughters
                  p->SetKS(11);

                  int nfp = plist->GetEntries();          // n. fragm. products
                  int ndp = decay_products->GetEntries(); // n. decay products

                  p->SetFirstChild ( nfp );          // decay products added at
                  p->SetLastChild  ( nfp + ndp -1 ); // the end of the fragm.rec.

                  //--  add decay products to the fragmentation record
                  TMCParticle * dp = 0;
                  TIter dpiter(decay_products);

                  while ( (dp = (TMCParticle *) dpiter.Next()) ) {

                     dp->SetParent(idecaying);
                     new ( (*plist)[plist->GetEntries()] ) TMCParticle(*dp);
                  }

                  //-- clean up decay products
                  decay_products->Delete();
                  delete decay_products;
               }

          } // particle is to be decayed
        } // KS < 10 : final state particle (as in PYTHIA LUJETS record)
     } // particles in fragmentation record
  } // force decay
}
//____________________________________________________________________________

