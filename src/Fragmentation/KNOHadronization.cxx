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

#include <TSystem.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TMCParticle6.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Decay/DecayModelI.h"
#include "Fragmentation/KNOHadronization.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils::print;

//____________________________________________________________________________
KNOHadronization::KNOHadronization() :
HadronizationModelBase("genie::KNOHadronization")
{
  fBaryonXFpdf  = 0;
  fBaryonPT2pdf = 0;
  fKNO          = 0;
}
//____________________________________________________________________________
KNOHadronization::KNOHadronization(string config) :
HadronizationModelBase("genie::KNOHadronization", config)
{
  fBaryonXFpdf  = 0;
  fBaryonPT2pdf = 0;
  fKNO          = 0;
}
//____________________________________________________________________________
KNOHadronization::~KNOHadronization()
{
  if (fBaryonXFpdf ) delete fBaryonXFpdf;
  if (fBaryonPT2pdf) delete fBaryonPT2pdf;
  if (fKNO         ) delete fKNO;
}
//____________________________________________________________________________
// HadronizationModelI interface implementation:
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

  if(!this->AssertValidity(interaction)) {
     LOG("KNOHad", pWARN) << "Returning a null particle list!";
     return 0;
  }

  double W = utils::kinematics::CalcW(interaction);
  LOG("KNOHad", pINFO) << "W = " << W << " GeV";

  //-- Select hadronic shower particles
  PDGCodeList * pdgcv = this->SelectParticles(interaction);

  //-- Decay the hadronic final state
  //   Two stratefies are considered (for N particles):
  //   1- N (>=2) particles get passed to the phase space decayer. This is the
  //      old NeuGEN strategy.
  //   2- The generated baryon P4 gets selected from from experimental xF and 
  //      pT^2 distributions and the remaining N-1 particles are passed to the
  //      phase space decayer, with P4 = P4(Sum_Hadronic) - P4(Baryon).
  //      For N=2, the meson P4 would be generated from the baryon P4 and 
  //      energy/momentrum conservation only. This is decay strategy adopted
  //      at the July-2006 hadronization model mini-workshop (C.Andreopoulos,
  //      H.Gallagher, T.Yang)

  TClonesArray * particle_list = 0;
  if(fUseBaryonXfPt2Param) particle_list = this->DecayMethod2(W,*pdgcv);
  else                     particle_list = this->DecayMethod1(W,*pdgcv);

  if(!particle_list) {
    LOG("KNOHad", pNOTICE) 
        << "Failed decaying a hadronic system @ W=" << W 
        	                    << "with  multiplicity=" << pdgcv->size();

    // clean-up and exit
    delete pdgcv;
    return 0;
  }

  //-- Handle unstable particle decays (if requested)
  this->HandleDecays(particle_list);

  //-- The container 'owns' its elements
  particle_list->SetOwner(true);

  return particle_list;
}
//____________________________________________________________________________
PDGCodeList * KNOHadronization::SelectParticles(
                                       const Interaction * interaction) const
{
  if(!this->AssertValidity(interaction)) {
     LOG("KNOHad", pWARN) << "Returning a null particle list!";
     return 0;
  }

  unsigned int min_mult = 2;
  unsigned int mult     = 0;
  PDGCodeList * pdgcv   = 0;

  double W = utils::kinematics::CalcW(interaction);

  //-- Get the charge that the hadron shower needs to have so as to
  //   conserve charge in the interaction
  int maxQ = this->HadronShowerCharge(interaction);
  LOG("KNOHad", pINFO) << "Hadron Shower Charge = " << maxQ;

   //-- Build the multiplicity probabilities for the input interaction
  LOG("KNOHad", pDEBUG) << "Building Multiplicity Probability distribution";
  LOG("KNOHad", pDEBUG) << *interaction;
  Option_t * opt = "+LowMultSuppr+Renormalize";
  TH1D * mprob = this->MultiplicityProb(interaction,opt);

  if(!mprob) {
    LOG("KNOHad", pERROR) << "Null multiplicity probability distribution!";
    return 0;
  }
  if(mprob->Integral("width")<=0) {
    LOG("KNOHad", pERROR) << "Empty multiplicity probability distribution!";
    delete mprob;
    return 0;
  }

  //----- FIND AN ALLOWED SOLUTION FOR THE HADRONIC FINAL STATE

  bool allowed_state=false;
  register unsigned int itry = 0;

  while(!allowed_state) 
  {
    itry++;

    //-- Go in error if a solution has not been found after many attempts
    if(itry>kMaxKNOHadSystIterations) {
       LOG("KNOHad", pERROR) 
         << "Couldn't select hadronic shower particles after: " 
                                                 << itry << " attempts!";
       delete mprob;
       return 0;
    }

    //-- Generate a hadronic multiplicity 
    mult = TMath::Nint( mprob->GetRandom() );

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
        delete mprob;
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
    //    bool permitted = (p4.Mag() > msum);
    bool permitted = (W > msum);

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

  delete mprob;

  return pdgcv;
}
//____________________________________________________________________________
TH1D * KNOHadronization::MultiplicityProb(
		        const Interaction * interaction, Option_t * opt) const
{
// Returns a multiplicity probability distribution for the input interaction.
// The input option (Default: "") can contain (combinations) of these strings:
//  - "+LowMultSuppr": applies NeuGEN Rijk factors suppresing the low multipl.
//    (1-pion and 2-pion) states as part of the DIS/RES joining scheme.
//  - "+Renormalize": renormalizes the probability distribution after applying
//    the NeuGEN scaling factors: Eg, when used as a hadronic multiplicity pdf 
//    the output hadronic multiplicity probability histogram needs to be re-
//    normalized. But, when this method is called from a DIS cross section 
//    algorithm using the integrated probability reduction as a cross section 
//    section reduction factor then the output histogram should not be re-
//    normalized after applying the scaling factors.

  if(!this->AssertValidity(interaction)) {
     LOG("KNOHad", pWARN) 
                << "Returning a null multiplicity probability distribution!";
     return 0;
  }

  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  // Compute the average charged hadron multiplicity as: <n> = a + b*ln(W^2)
  // Calculate avergage hadron multiplicity (= 1.5 x charged hadron mult.)

  double W     = utils::kinematics::CalcW(interaction);
  double avnch = this->AverageChMult(nu_pdg, nuc_pdg, W);
  double avn   = 1.5*avnch;

  SLOG("KNOHad", pINFO) 
             << "Average hadronic multiplicity (W=" << W << ") = " << avn;

  // Find the max possible multiplicity as W = Mneutron + (maxmult-1)*Mpion
  double maxmult = this->MaxMult(interaction);;

  // If required force the NeuGEN maximum multiplicity limit (10)
  // Note: use for NEUGEN/GENIE comparisons, not physics MC production
  if(fForceNeuGenLimit & maxmult>10) maxmult=10;

  // Set maximum multiplicity so that it does not exceed the max number of
  // particles accepted by the ROOT phase space decayer (18)
  // Change this if ROOT authors remove the TGenPhaseSpace limitation.
  if(maxmult>18) maxmult=18;

  SLOG("KNOHad", pDEBUG) << "Computed maximum multiplicity = " << maxmult;

  // Create multiplicity probability histogram
  TH1D * mult_prob = this->CreateMultProbHist(maxmult);

  // Compute the multiplicity probabilities values up to the bin corresponding 
  // to the computed maximum multiplicity

  if(maxmult>2) {
    int nbins = mult_prob->FindBin(maxmult);

    for(int i = 1; i <= nbins; i++) {
       // KNO distribution is <n>*P(n) vs n/<n>
       double n    = mult_prob->GetBinCenter(i);  // bin centre
       double z    = n/avn;                       // z=n/<n>
       double avnP = this->KNO(nu_pdg,nuc_pdg,z); // <n>*P(n)
       double P    = avnP / avn;                  // P(n)

       SLOG("KNOHad", pDEBUG)
          << "n = " << n << " (n/<n> = " << z
                            << ", <n>*P = " << avnP << ") => P = " << P;
       mult_prob->Fill(n,P);
    }
  } else {
       mult_prob->Fill(2,1.);
  }

  double integral = mult_prob->Integral("width");
  if(integral>0) {
    // Normalize the probability distribution
    mult_prob->Scale(1.0/integral);
  } else {
    SLOG("KNOHad", pWARN) << "probability distribution integral = 0";
    return mult_prob;
  }

  string option(opt);

  bool apply_neugen_Rijk = option.find("+LowMultSuppr") != string::npos;
  bool renormalize       = option.find("+Renormalize")  != string::npos;

  // Apply the NeuGEN probability scaling factors -if requested-
  if(apply_neugen_Rijk) {
    SLOG("KNOHad", pINFO) << "Applying NeuGEN scaling factors";
     // Only do so for W<Wcut
     if(W<fWcut) {
       this->ApplyRijk(interaction, renormalize, mult_prob);
     } else {
        SLOG("KNOHad", pDEBUG)  
              << "W = " << W << " < Wcut = " << fWcut 
                                << " - Will not apply scaling factors";
     }//<wcut?
  }//apply?

  return mult_prob;
}
//____________________________________________________________________________
double KNOHadronization::Weight(void) const
{
  return fWeight;
}
//____________________________________________________________________________
// methods overloading the default Algorithm interface implementation:
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
// private methods:
//____________________________________________________________________________
void KNOHadronization::LoadConfig(void)
{
// Read configuration options or set defaults

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Delete the KNO spline from previous configuration of this instance
  if(fKNO) delete fKNO;

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

  // Decay unstable particles now or leave it for later? Which decayer to use?
  fDecayer = 0;
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

  // load legacy KNO spline
  fUseLegacyKNOSpline = fConfig->GetBoolDef("use-legacy-kno-spline", false);

  if(fUseLegacyKNOSpline) {
     assert(gSystem->Getenv("GENIE"));
     string basedir    = gSystem->Getenv("GENIE");
     string defknodata = basedir + "/data/kno/KNO.dat";
     string knodata    = fConfig->GetStringDef("kno-data",defknodata);

     LOG("KNOHad", pNOTICE) << "Loading KNO data from: " << knodata;

     fKNO = new Spline(knodata);
  }

  // Load parameters determining the average multiplicity
  fAvp  = fConfig->GetDoubleDef("alpha-vp",  gc->GetDouble("KNO-Alpha-vp") );
  fAvn  = fConfig->GetDoubleDef("alpha-vn",  gc->GetDouble("KNO-Alpha-vn") ); 
  fAvbp = fConfig->GetDoubleDef("alpha-vbp", gc->GetDouble("KNO-Alpha-vbp")); 
  fAvbn = fConfig->GetDoubleDef("alpha-vbn", gc->GetDouble("KNO-Alpha-vbn"));
  fB    = fConfig->GetDoubleDef("beta",      gc->GetDouble("KNO-Beta")     );

  // Load the Levy function parameter
  fCvp  = fConfig->GetDoubleDef("levy-c-vp",  gc->GetDouble("KNO-LevyC-vp") );
  fCvn  = fConfig->GetDoubleDef("levy-c-vn",  gc->GetDouble("KNO-LevyC-vn") ); 
  fCvbp = fConfig->GetDoubleDef("levy-c-vbp", gc->GetDouble("KNO-LevyC-vbp")); 
  fCvbn = fConfig->GetDoubleDef("levy-c-vbn", gc->GetDouble("KNO-LevyC-vbn"));

  // Force NEUGEN upper limit in hadronic multiplicity (to be used only
  // NEUGEN/GENIE comparisons)
  fForceNeuGenLimit = fConfig->GetBoolDef("force-neugen-mult-limit", false);

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));

  // Load NEUGEN multiplicity probability scaling parameters Rijk
  fRvpCCm2  = fConfig->GetDoubleDef(
                      "R-vp-CC-m2", gc->GetDouble("DIS-HMultWgt-vp-CC-m2"));
  fRvpCCm3  = fConfig->GetDoubleDef(
                      "R-vp-CC-m3", gc->GetDouble("DIS-HMultWgt-vp-CC-m3"));
  fRvpNCm2  = fConfig->GetDoubleDef(
                      "R-vp-NC-m2", gc->GetDouble("DIS-HMultWgt-vp-NC-m2"));
  fRvpNCm3  = fConfig->GetDoubleDef(
                      "R-vp-NC-m3", gc->GetDouble("DIS-HMultWgt-vp-NC-m3"));
  fRvnCCm2  = fConfig->GetDoubleDef(
                      "R-vn-CC-m2", gc->GetDouble("DIS-HMultWgt-vn-CC-m2"));
  fRvnCCm3  = fConfig->GetDoubleDef(
                      "R-vn-CC-m3", gc->GetDouble("DIS-HMultWgt-vn-CC-m3"));
  fRvnNCm2  = fConfig->GetDoubleDef(
                      "R-vn-NC-m2", gc->GetDouble("DIS-HMultWgt-vn-NC-m2"));
  fRvnNCm3  = fConfig->GetDoubleDef(
                      "R-vn-NC-m3", gc->GetDouble("DIS-HMultWgt-vn-NC-m3"));
  fRvbpCCm2 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m2",gc->GetDouble("DIS-HMultWgt-vbp-CC-m2"));
  fRvbpCCm3 = fConfig->GetDoubleDef(
                     "R-vbp-CC-m3",gc->GetDouble("DIS-HMultWgt-vbp-CC-m3"));
  fRvbpNCm2 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m2",gc->GetDouble("DIS-HMultWgt-vbp-NC-m2"));
  fRvbpNCm3 = fConfig->GetDoubleDef(
                     "R-vbp-NC-m3",gc->GetDouble("DIS-HMultWgt-vbp-NC-m3"));
  fRvbnCCm2 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m2",gc->GetDouble("DIS-HMultWgt-vbn-CC-m2"));
  fRvbnCCm3 = fConfig->GetDoubleDef(
                     "R-vbn-CC-m3",gc->GetDouble("DIS-HMultWgt-vbn-CC-m3"));
  fRvbnNCm2 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m2",gc->GetDouble("DIS-HMultWgt-vbn-NC-m2"));
  fRvbnNCm3 = fConfig->GetDoubleDef(
                     "R-vbn-NC-m3",gc->GetDouble("DIS-HMultWgt-vbn-NC-m3"));

  LOG("KNOHad", pINFO) << *fConfig;
}
//____________________________________________________________________________
double KNOHadronization::KNO(int nu_pdg, int nuc_pdg, double z) const
{
// Computes <n>P(n) for the input reduced multiplicity z=n/<n>

  if(fUseLegacyKNOSpline) {
     bool inrange = z > fKNO->XMin() && z < fKNO->XMax();
     return (inrange) ? fKNO->Evaluate(z) : 0.;
  }
  assert( pdg::IsNeutrino(nu_pdg) || pdg::IsAntiNeutrino(nu_pdg) );
  assert( pdg::IsNeutronOrProton(nuc_pdg) );

  double c=0; // Levy function parameter

  if      ( pdg::IsProton (nuc_pdg) && pdg::IsNeutrino    (nu_pdg) ) c=fCvp;
  else if ( pdg::IsNeutron(nuc_pdg) && pdg::IsNeutrino    (nu_pdg) ) c=fCvn;
  else if ( pdg::IsProton (nuc_pdg) && pdg::IsAntiNeutrino(nu_pdg) ) c=fCvbp;
  else if ( pdg::IsNeutron(nuc_pdg) && pdg::IsAntiNeutrino(nu_pdg) ) c=fCvbn;
  else return 0;

  double x   = c*z+1;
  double kno = 2*TMath::Exp(-c)*TMath::Power(c,x)/TMath::Gamma(x);

  return kno;
}
//____________________________________________________________________________
double KNOHadronization::AverageChMult(int nu_pdg,int nuc_pdg, double W) const
{
// computes the average charged multiplicity
//
  assert( pdg::IsNeutrino(nu_pdg) || pdg::IsAntiNeutrino(nu_pdg) );
  assert( pdg::IsNeutronOrProton(nuc_pdg) );

  double a=0, b=fB;

  if      ( pdg::IsProton (nuc_pdg) && pdg::IsNeutrino    (nu_pdg) ) a=fAvp;
  else if ( pdg::IsNeutron(nuc_pdg) && pdg::IsNeutrino    (nu_pdg) ) a=fAvn;
  else if ( pdg::IsProton (nuc_pdg) && pdg::IsAntiNeutrino(nu_pdg) ) a=fAvbp;
  else if ( pdg::IsNeutron(nuc_pdg) && pdg::IsAntiNeutrino(nu_pdg) ) a=fAvbn;
  else return 0;

  double av_nch = a + b * 2*TMath::Log(W);
  return av_nch;        
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
TClonesArray * KNOHadronization::DecayMethod1(
                                     double W, const PDGCodeList & pdgv) const
{
// Simple phase space decay including all generated particles.
// The old NeuGEN decay strategy.

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
                                     double W, const PDGCodeList & pdgv) const
{
// Generate the baryon based on experimental pT^2 and xF distributions
// Then pass the remaining system of N-1 particles to a phase space decayer.
// The strategy adopted at the July-2006 hadronization model mini-workshop.

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
  bool allowdup = true;
  PDGCodeList pdgv_stripped(pdgv.size()-1, allowdup);
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
                                     double W, const PDGCodeList & pdgv) const
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
        	                   const PDGCodeList & pdgv, int offset) const
{
// General method decaying the input particle system 'pdgv' with available 4-p
// given by 'pd'. The decayed system is used to populate the input TMCParticle 
// array starting from the slot 'offset'.
//
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
double KNOHadronization::ReWeightPt2(const PDGCodeList & pdgcv) const
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
PDGCodeList * KNOHadronization::GenerateFSHadronCodes(
                                   int multiplicity, int maxQ, double W) const
{
// Selection of fragments (identical as in NeuGEN).

  // PDG Library
  PDGLibrary * pdg = PDGLibrary::Instance();

  // vector to add final state hadron PDG codes
  bool allowdup=true;
  PDGCodeList * pdgc = new PDGCodeList(allowdup);
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
// Selection of main target fragment (identical as in NeuGEN).
// Assign baryon as p or n. Force it for ++ and - I=3/2 at mult. = 2

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
bool KNOHadronization::AssertValidity(const Interaction * interaction) const
{
  if(interaction->GetExclusiveTag().IsCharmEvent()) {
     LOG("KNOHad", pWARN) << "Can't hadronize charm events";
     return false;
  }

  double W = utils::kinematics::CalcW(interaction);
  if(W < this->Wmin()) {
     LOG("KNOHad", pWARN) << "Low invariant mass, W = " << W << " GeV!!";
     return false;
  }
  return true;
}
//____________________________________________________________________________


