//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    The KNO hadronization model.

          This hadronization scheme is similar to the one originally used
          in NeuGEN by G.Barr, G.F.Pearce, H.Gallagher. \n

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TMCParticle6.h>
#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Decay/DecayModelI.h"
#include "Fragmentation/KNOHadronization.h"
#include "Fragmentation/MultiplicityProb.h"
#include "Fragmentation/MultiplicityProbModelI.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
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
  //----- Create a multiplicity probability object & attach the requested model

  LOG("KNOHad", pDEBUG) << "Creating MultiplicityProb & attaching model";

  MultiplicityProb mult_prob;
  mult_prob.AttachModel(fMultProbModel);

  //----- Build the multiplicity probabilities for the input interaction

  LOG("KNOHad", pDEBUG) << "Building Multiplicity Probability distribution";
  LOG("KNOHad", pDEBUG) << *interaction;

  mult_prob.BuildProbDistribution(interaction);

  //----- Get a random multiplicity based of their probability distribution
  unsigned int min_mult = 2;
  unsigned int max_mult = 20;
  unsigned int mult = mult_prob.RandomMultiplicity(min_mult, max_mult);

  //----- Get the charge that the hadron shower needs to have so as to
  //      conserve charge in the interaction
  int maxQ = HadronShowerCharge(interaction);

  //----- force generated multiplicty to be consistent with the
  //      hadron shower charge needed to be generated
  //      eg. if mult = 1 but maxQ = 2, force mult = 2
  mult = TMath::Max(mult, (unsigned int) TMath::Abs(maxQ));

  LOG("KNOHad", pINFO) << "Hadron multiplicity  = " << mult;
  LOG("KNOHad", pINFO) << "Hadron Shower Charge = " << maxQ;

  if(mult < min_mult || mult > max_mult) {
     LOG("KNOHad", pERROR) 
                      << "Multiplicity out of bounds. Return NULL list";
     return 0;
  }

  //----- Initial particle 4-momentum
  //      (All energy = final state invariant mass)
  double W = interaction->GetKinematics().W();
  LOG("KNOHad", pINFO) << "W = " << W << " GeV";

  assert(W > kNucleonMass+kPionMass);

  TLorentzVector p4(0,0,0,W);

  //----- Determine what kind of particles we have in the final state

  // This section was translated / adapted from NeuGEN.
  // Original author: H.Gallagher
  vector<int> * pdgc = GenerateFSHadronCodes(mult, maxQ, W);

  LOG("KNOHad", pDEBUG)
             << "Generated multiplicity: " << pdgc->size();

  // muliplicity might have been forced to smaller value if the invariant
  // mass of the hadronic system was not sufficient
  mult = pdgc->size(); // update for potential change

  //----- DECAY INCLUSIVES

  // Determine the kinematics of the final state particles using
  // ROOT's Phase Space Generator

  LOG("KNOHad", pDEBUG) << "Generating phase space";

  vector<int>::const_iterator pdg_iter;
  double mass[pdgc->size()];
  int i = 0;
  for(pdg_iter = pdgc->begin(); pdg_iter != pdgc->end(); ++pdg_iter)
          mass[i++] = PDGLibrary::Instance()->Find(*pdg_iter)->Mass();

  bool permitted = fPhaseSpaceGenerator.SetDecay(p4, mult, mass);

  if(!permitted) {
     LOG("KNOHad", pERROR) << "*** Decay forbidden by kinematics! ***";
     LOG("KNOHad", pERROR) << "This should not have happened";
     LOG("KNOHad", pERROR) << "Initial State: ";
     LOG("KNOHad", pERROR)
            << "P4 = (" << p4.Px() << ", " << p4.Py()
            <<     ", " << p4.Pz() << ", " << p4.Energy() << ")";
     LOG("KNOHad", pERROR) << "Final State: ";
     for(unsigned int i = 0; i < pdgc->size(); i++) {
       LOG("KNOHad", pERROR)
        << " -PDGC = " << (*pdgc)[i] << ", m = " << mass[i] << " GeV";
     }
     return 0;
  }

  //-- generate kinematics in the Center-of-Mass (CM) frame
  fPhaseSpaceGenerator.Generate();

  //-- insert final state products into a TClonesArray of TMCParticles
  //   and return it
  LOG("KNOHad", pDEBUG)
       << "Creating output particle list TClonesArray of TMCParticles";

  TClonesArray * particle_list = new TClonesArray("TMCParticle", mult);

  for(unsigned int i = 0; i < pdgc->size(); i++) {

     //-- get the 4-momentum of the i-th final state particle
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(i);

     //-- push TMCParticle in the particle list (TClonesArray)
     LOG("KNOHad", pDEBUG)
               << "Adding final state particle PDGC = " << (*pdgc)[i]
               << " with mass = " << mass[i] << " GeV";

     new ( (*particle_list)[i] ) TMCParticle
        (
           1,               /* KS Code                          */
           (*pdgc)[i],      /* PDG Code                         */
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

  delete pdgc;

  return particle_list;
}
//____________________________________________________________________________
void KNOHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void KNOHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void KNOHadronization::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed

  fMultProbModel = 0;
  fDecayer       = 0;

  fMultProbModel = dynamic_cast<const MultiplicityProbModelI *> (
    this->SubAlg("multiplicity-prob-alg-name", "multiplicity-prob-param-set"));
  assert(fMultProbModel);

  if(fForceDecays) {
      fDecayer = dynamic_cast<const DecayModelI *> (
                       this->SubAlg("decayer-alg-name", "decayer-param-set"));
      assert(fDecayer);
  }
}
//____________________________________________________________________________
void KNOHadronization::LoadConfigData(void)
{
  fForceDecays = fConfig->GetBoolDef("force-decays", false);

  // Probability for producing hadron each pairs
  fPpi0 = fConfig->GetDoubleDef("prob-fs-pi0-pair",       0.30); // pi0 pi0
  fPpic = fConfig->GetDoubleDef("prob-fs-piplus-piminus", 0.60); // pi+ pi-
  fPKc  = fConfig->GetDoubleDef("prob-fs-Kplus-Kminus",   0.05); // K+  K-
  fPK0  = fConfig->GetDoubleDef("prob-fs-K0-K0bar",       0.05); // K0  K0bar
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
        LOG("KNOHad", pDEBUG)
                  << "\n Need more negative charge ---> Adding a pi-";
        pdgc->push_back( kPdgPiMinus );

        // update n-of-hadrons to add, avail. shower charge & invariant mass
        maxQ += 1;
        hadrons_to_add--;

        W -= pdg->Find(kPdgPiMinus)->Mass();

     } else if (maxQ > 0) {
        // Need more positive charge
        LOG("KNOHad", pDEBUG)
                  << "\n Need more positive charge ---> Adding a pi+";
        pdgc->push_back( kPdgPiPlus );

        // update n-of-hadrons to add, avail. shower charge & invariant mass
        maxQ -= 1;
        hadrons_to_add--;

        W -= pdg->Find(kPdgPiPlus)->Mass();
     }
  }

  //-- Add remaining neutrals or pairs up to the generated multiplicity
  if(maxQ == 0) {

     LOG("KNOHad", pDEBUG) <<
      "\n Final state has correct Q - Now adding only neutrals or +- pairs";

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
               << "\n Odd number of hadrons left to add ---> Adding a pi0";
        pdgc->push_back( kPdgPi0 );

        // update n-of-hadrons to add & available invariant mass
        hadrons_to_add--;
        W -= pdg->Find(kPdgPi0)->Mass();
     }

     // Now add pairs (pi0 pi0 / pi+ pi- / K+ K- / K0 K0bar)

     assert( hadrons_to_add % 2 == 0 ); // even number

     RandomGen * rnd = RandomGen::Instance();

     while(hadrons_to_add > 0 && W >= M2pi0) {

         double x = rnd->Random1().Rndm();
         LOG("KNOHad", pDEBUG) << "rndm = " << x;

         if (x >= 0 && x < fPpi0) {
            //-------------------------------------------------------------
            // Add a pi0-pair
            LOG("KNOHad", pDEBUG) << "\n ---> Adding a pi0pi0 pair";

            pdgc->push_back( kPdgPi0 );
            pdgc->push_back( kPdgPi0 );

            hadrons_to_add -= 2; // update the number of hadrons to add
            W -= M2pi0; // update the available invariant mass
            //-------------------------------------------------------------

         } else if (x < fPpi0 + fPpic) {
            //-------------------------------------------------------------
            // Add a piplus - piminus pair if there is enough W
            if(W >= M2pic) {

                LOG("KNOHad", pDEBUG) << "\n ---> Adding a pi+pi- pair";

                pdgc->push_back( kPdgPiPlus  );
                pdgc->push_back( kPdgPiMinus );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2pic; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) <<
                  "\n Not enough mass for a pi+pi-: trying something else";
            }
            //-------------------------------------------------------------

         } else if (x < fPpi0 + fPpic + fPKc) {
            //-------------------------------------------------------------
            // Add a Kplus - Kminus pair if there is enough W
            if(W >= M2Kc) {

                LOG("KNOHad", pDEBUG) << "\n ---> Adding a K+K- pair";

                pdgc->push_back( kPdgKPlus  );
                pdgc->push_back( kPdgKMinus );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2Kc; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) <<
                    "\n Not enough mass for a K+K-: trying something else";
            }
            //-------------------------------------------------------------

         } else if (x <= fPpi0 + fPpic + fPKc + fPK0) {
            //-------------------------------------------------------------
            // Add a K0 - K0bar pair if there is enough W
            if( W >= M2K0 ) {
                LOG("KNOHad", pDEBUG) << "\n ---> Adding a K0 K0bar pair";

                pdgc->push_back( kPdgK0 );
                pdgc->push_back( kPdgK0 );

                hadrons_to_add -= 2; // update the number of hadrons to add
                W -= M2K0; // update the available invariant mass
            } else {
                LOG("KNOHad", pDEBUG) <<
                 "\n Not enough mass for a K0 K0bar: trying something else";
            }
            //-------------------------------------------------------------

         } else {
            LOG("KNOHad", pERROR)
                 << "FS Hadron Assignment Probabilities do not add up to 1!!";
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
  double x = rnd->Random1().Rndm();

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

  LOG("KNOHad", pDEBUG) <<
        "\n ---> Adding a " << (( pdgc == kPdgProton ) ? "proton" : "neutron");

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

