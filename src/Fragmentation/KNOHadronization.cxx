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
#include <TF1.h>

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
  fBaryonXFpdf  = 0;
  fBaryonPT2pdf = 0;
}
//____________________________________________________________________________
KNOHadronization::KNOHadronization(string config) :
HadronizationModelI("genie::KNOHadronization", config)
{
  fBaryonXFpdf  = 0;
  fBaryonPT2pdf = 0;
}
//____________________________________________________________________________
KNOHadronization::~KNOHadronization()
{
  if (fBaryonXFpdf ) delete fBaryonXFpdf;
  if (fBaryonPT2pdf) delete fBaryonPT2pdf;
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

  unsigned int min_mult = 2;
  unsigned int mult     = 0;
  vector<int> * pdgcv   = 0;

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

    // is it an allowed decay?
    double msum=0;
    vector<int>::const_iterator pdg_iter;
    for(pdg_iter = pdgcv->begin(); pdg_iter != pdgcv->end(); ++pdg_iter) {
      int pdgc = *pdg_iter;
      double m = PDGLibrary::Instance()->Find(pdgc)->Mass();

      msum += m;
      LOG("KNOHad", pDEBUG) << "- PDGC=" << pdgc << ", m=" << m << " GeV";
    }
    bool permitted = (p4.Mag() > msum);

    if(!permitted) {
       LOG("KNOHad", pWARN) << "*** Decay forbidden by kinematics! ***";
       LOG("KNOHad", pWARN) << "sum{mass} = " << msum << ", W = " << W;
       LOG("KNOHad", pWARN) << "Discarding hadronic system & re-trying!";
       delete pdgcv;
       allowed_state = false;
       continue;
    }

    allowed_state = true;
    LOG("KNOHad", pNOTICE) 
             << "Found an allowed hadronic state @ W=" << W 
                                              << " multiplicity=" << mult;
  } // attempts


  //----- DECAY HADRONIC FINAL STATE

  //-- Decide the hadronic system decay strategy
  //   Options considered for N particles:
  //   1- N (>=2) particles get passed to the phase space decayer
  //   2- The generated baryon P4 gets selected from from experimental xF and pT^2 
  //      distributions and the remaining N-1 particles are passed to the phase
  //      space decayer, with P4 = P4(Sum_Hadronic) - P4(Baryon)
  //      For N=2, the meson P4 would be generated from the baryon P4 and energy/
  //      momentrum conservation only.

  TClonesArray * particle_list = 0;
  if(fUseBaryonXfPt2Param) particle_list = this->DecayMethod2(W,*pdgcv);
  else                     particle_list = this->DecayMethod1(W,*pdgcv);

  if(!particle_list) {
    LOG("KNOHad", pNOTICE) 
        << "Failed decaying a hadronic system @ W=" << W 
                                             << "with  multiplicity=" << mult;

    // clean-up and exit
    delete pdgcv;
    return 0;
  }

  // handle unstable particle decays (if requested)
  this->HandleDecays(particle_list);

  // the container 'owns' its elements
  particle_list->SetOwner(true);

  return particle_list;
}
//____________________________________________________________________________
TClonesArray * KNOHadronization::DecayMethod1(
                                     double W, const vector<int> & pdgv) const
{
// Simple phase space decay including all generated particles
//

  LOG("KNOHad", pINFO) << "** Using Hadronic System Decay method 1";

  TLorentzVector p4had(0,0,0,W);
  TClonesArray * plist = new TClonesArray("TMCParticle", pdgv.size());

  bool ok = this->PhaseSpaceDecay(*plist, p4had, pdgv); // do the decay

  // clean-up and return NULL
  if(!ok) {
     plist->Delete();
     delete plist;
     return 0;
  }
  return plist;
}
//____________________________________________________________________________
TClonesArray * KNOHadronization::DecayMethod2(
                                     double W, const vector<int> & pdgv) const
{
// Generate the baryon based on experimental pT^2 and xF distributions
// Then pass the remaining system of N-1 particles to a phase space decayer

  LOG("KNOHad", pINFO) << "** Using Hadronic System Decay method 2";

  // If only 2 particles are input then don't call the phase space decayer
  if(pdgv.size() == 2) return this->DecayBackToBack(W,pdgv);

  // Now handle the more general case:

  // Take the baryon
  int    baryon = pdgv[0]; 
  double MN     = PDGLibrary::Instance()->Find(baryon)->Mass();
  double MN2    = TMath::Power(MN, 2);

  assert(pdg::IsNeutronOrProton(baryon));

  // Strip the PDG list from the baryon
  vector<int> pdgv_stripped(pdgv.size()-1);
  for(unsigned int i=1; i<pdgv.size(); i++) pdgv_stripped[i-1] = pdgv[i];
  
  // Get the sum of all masses for the particles in the stripped list

  double mass_sum = 0;
  vector<int>::const_iterator pdg_iter = pdgv_stripped.begin();

  for( ; pdg_iter != pdgv_stripped.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    mass_sum += PDGLibrary::Instance()->Find(pdgc)->Mass();
  }

  // generate the N 4-p independently

  LOG("KNOHad", pINFO) << "Generating p4 for baryon with pdg= " << baryon;

  RandomGen * rnd = RandomGen::Instance();
  TLorentzVector p4had(0,0,0,W);
  TLorentzVector p4N  (0,0,0,0);
  TLorentzVector p4d;

  bool allowed = false;

  while(!allowed) {

    //-- generate baryon xF and pT2
    double xf  = fBaryonXFpdf->GetRandom();               
    double pt2 = fBaryonPT2pdf->GetRandom();               

    //-- generate baryon px,py,pz
    double pt  = TMath::Sqrt(pt2);            
    double phi = (2*kPi) * rnd->RndHadro().Rndm();
    double px  = pt * TMath::Cos(phi);
    double py  = pt * TMath::Sin(phi);
    double pz  = xf*W/2;
    double p2  = TMath::Power(pz,2) + pt2;
    double E   = TMath::Sqrt(p2+MN2);

    p4N.SetPxPyPzE(px,py,pz,E);

    LOG("KNOHad", pDEBUG) << "Trying nucleon xF= "<< xf<< ", pT2= "<< pt2;

    //-- check whether there is phase space for the remnant N-1 system
   
    p4d = p4had-p4N; // 4-momentum vector for phase space decayer

    double Mav = p4d.Mag();
    allowed = (Mav > mass_sum);

    if(allowed) {
      LOG("KNOHad", pINFO) 
        << "Generated baryon with P4 = " << utils::print::P4AsString(&p4N);
      LOG("KNOHad", pINFO) 
        << "Remaining available mass = " << Mav 
                                     << ", particle masses = " << mass_sum; 
    }
  }

  // Create the particle list
  TClonesArray * plist = new TClonesArray("TMCParticle", pdgv.size());

  // Insert the baryon
  new ((*plist)[0]) TMCParticle(
    1,baryon,-1,-1,-1, p4N.Px(),p4N.Py(),p4N.Pz(),p4N.Energy(),MN, 0,0,0,0,0);

  // Do a phase space decay for the N-1 particles and add them to the list
  bool ok = this->PhaseSpaceDecay(*plist, p4d, pdgv_stripped, 1);

  // clean-up and return NULL
  if(!ok) {
     plist->Delete();
     delete plist;
     return 0;
  }
  return plist;
}
//____________________________________________________________________________
TClonesArray * KNOHadronization::DecayBackToBack(
                                     double W, const vector<int> & pdgv) const
{
// Handles a special case (only two particles) of the 2nd decay method 
//
  LOG("KNOHad", pINFO) << "Generating two particles back-to-back";

  assert(pdgv.size()==2);

  // Get the particle PDG codes and masses
  int    baryon = pdgv[0]; 
  int    meson  = pdgv[1]; 
  double MN     = PDGLibrary::Instance()->Find(baryon)->Mass();
  double Mpi    = PDGLibrary::Instance()->Find(meson)->Mass();
  double MN2    = TMath::Power(MN, 2);
 
  assert(pdg::IsNeutronOrProton(baryon));

  // Create the particle list

  TClonesArray * plist = new TClonesArray("TMCParticle", pdgv.size());

  // Do the decay

  RandomGen * rnd = RandomGen::Instance();
  bool accepted = false;
 
  double Ro = 1.2 * fBaryonPT2pdf->GetMaximum(0,1);

  while(!accepted) {

    //-- generate baryon xF and pT2
    double xfN = fBaryonXFpdf->GetRandom();

    //-- compute pL (=pz)
    double pz  = xfN*W/2;
    double pz2 = TMath::Power(pz,2);

    //-- generate EN in (0,W)
    double EN  = W * rnd->RndHadro().Rndm();
    double EN2 = TMath::Power(EN,2);

    //-- is it allowed? (pT2>=0)
    double pT2  = EN2 - MN2 - pz2; 

    if(pT2<=0) continue; // reject current selection

    //-- use the rejection method to decide whether to accept the pT2 value

    double Rrnd = Ro * rnd->RndHadro().Rndm();
    double Rpdf = fBaryonPT2pdf->Eval(pT2);               

    accepted = (Rrnd < Rpdf);

    if(accepted) {
      // Finish up & insert the particles at the particle list

      double pT  = TMath::Sqrt(pT2);            
      double phi = (2*kPi) * rnd->RndHadro().Rndm();
      double px  = pT * TMath::Cos(phi);
      double py  = pT * TMath::Sin(phi);
      double Epi = W-EN;

      TClonesArray & particle_list = *plist;

      new (particle_list[0]) 
                TMCParticle (1,baryon,-1,-1,-1, px, py, pz,EN, MN, 0,0,0,0,0);
      new (particle_list[1]) 
                TMCParticle (1,meson, -1,-1,-1,-px,-py,-pz,Epi,Mpi,0,0,0,0,0);
    }
  }
  return plist;
}
//____________________________________________________________________________
bool KNOHadronization::PhaseSpaceDecay(
	   TClonesArray & plist, TLorentzVector & pd, 
       	                          const vector<int> & pdgv, int offset) const
{
  LOG("KNOHad", pINFO) << "*** Performing a Phase Space Decay";

  assert ( offset      >= 0);
  assert ( pdgv.size() >  1);

  // Get the decay product masses

  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  LOG("KNOHad", pINFO)  
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;
  LOG("KNOHad", pINFO) 
    << "Decaying system p4 = " << utils::print::P4AsString(&pd);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(pd, pdgv.size(), mass);
  if(!permitted) {
     LOG("KNOHad", pERROR) 
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(&pd);

     // clean-up and return
     delete [] mass;
     return false;
  }

  // Get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int i=0; i<200; i++) {
     double w = fPhaseSpaceGenerator.Generate();   
     w *= this->ReWeightPt2(pdgv);
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);

  LOG("KNOHad", pNOTICE) 
     << "Max phase space gen. weight @ current hadronic system: " << wmax;

  // Generate a weighted or unweighted decay

  RandomGen * rnd = RandomGen::Instance();

  if(fGenerateWeighted) 
  {
    // *** generating weighted decays ***
    double w = fPhaseSpaceGenerator.Generate();   
    w *= this->ReWeightPt2(pdgv); 
    fWeight *= TMath::Max(w/wmax, 1.);
  }
  else 
  {
    // *** generating un-weighted decays ***
     wmax *= 1.2;
     bool accept_decay=false;
     unsigned int itry=0;

     while(!accept_decay) 
     {
       itry++;

       if(itry>kMaxUnweightDecayIterations) {
         // report, clean-up and return
         LOG("KNOHad", pWARN) 
             << "Couldn't generate an unweighted phase space decay after " 
             << itry << " attempts";
         delete [] mass;
         return false;
       }

       double w  = fPhaseSpaceGenerator.Generate();   
       w *= this->ReWeightPt2(pdgv); 
       double gw = wmax * rnd->RndHadro().Rndm();

       LOG("KNOHad", pINFO) << "Decay weight = " << w << " / R = " << gw;

       accept_decay = (gw<=w);
     }
  }

  // Insert final state products into a TClonesArray of TMCParticles

  i=0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {

     //-- current PDG code
     int pdgc = *pdg_iter;

     //-- get the 4-momentum of the i-th final state particle
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(i);

     new ( plist[offset+i] ) TMCParticle(
           1,               /* KS Code                          */
           pdgc,            /* PDG Code                         */
          -1,               /* parent particle                  */
          -1,               /* first child particle             */
          -1,               /* last child particle              */
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
     i++;
  }

  // Clean-up
  delete [] mass;

  return true;
}
//____________________________________________________________________________
double KNOHadronization::ReWeightPt2(const vector<int> & pdgcv) const
{
// Phase Space Decay re-weighting to reproduce exp(-pT2/<pT2>) pion pT2 
// distributions.
// See: A.B.Clegg, A.Donnachie, A Description of Jet Structure by pT-limited
// Phase Space.

  if(!fReWeightDecays) return 1.;

  double w = 1;

  for(unsigned int i = 0; i < pdgcv.size(); i++) {

     int pdgc = pdgcv[i];

     if(pdgc!=kPdgPiPlus&&pdgc!=kPdgPiMinus) continue;

     TLorentzVector * p4 = fPhaseSpaceGenerator.GetDecay(i); 
     double pt2 = TMath::Power(p4->Px(),2) + TMath::Power(p4->Py(),2);
     double wi  = TMath::Exp(-fPhSpRwA*TMath::Sqrt(pt2));
     //double wi = (9.41 * TMath::Landau(pt2,0.24,0.12));

     w *= wi;
  }
  return w;
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

  // Generate the baryon xF and pT^2 using experimental data as PDFs? 
  // In this case, only the N-1 other particles would be fed into the phase
  // space decayer. This seems to improve hadronic system features such as 
  // bkw/fwd xF hemisphere average multiplicities.
  // Note: not in the legacy KNO model (NeuGEN). Switch this feature off for 
  // comparisons or for reproducing old simulations.
  fUseBaryonXfPt2Param = fConfig->GetBoolDef("use-baryon-xF-pT2-parm", true);

  // Reweight the phase space decayer events to reproduce the experimentally
  // measured pT^2 distributions.
  // Note: not in the legacy KNO model (NeuGEN). Switch this feature off for 
  // comparisons or for reproducing old simulations.
  fReWeightDecays = fConfig->GetBoolDef("reweight-phase-space-decays", true);

  // Generated weighted or un-weighted hadronic systems?
  fGenerateWeighted = fConfig->GetBoolDef("generate-weighted", false);

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

  // Baryon pT^2 and xF parameterizations used as PDFs

  if (fBaryonXFpdf ) delete fBaryonXFpdf;
  if (fBaryonPT2pdf) delete fBaryonPT2pdf;

  fBaryonXFpdf  = new TF1("fBaryonXFpdf",
                   "0.049889-0.116769*x-0.033949*x*x+0.146209*x*x*x",-1,0.5);  
  fBaryonPT2pdf = new TF1("fBaryonPT2pdf", "exp(-0.213362-6.62464*x)",0,0.6);  

  // Parameter for phase space re-weighting. See ReWeightPt2()

  fPhSpRwA = fConfig->GetDoubleDef(
           "phase-space-reweighting-param", gc->GetDouble("KNO-PhSpRwParm")); 

  LOG("KNOHad", pINFO) << *fConfig;
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

  int baryon_code = this->GenerateBaryonPdgCode(multiplicity, maxQ);
  pdgc->push_back(baryon_code);

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

