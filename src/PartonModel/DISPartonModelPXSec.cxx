//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 05, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <TH1D.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/DISStructureFuncModelI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Fragmentation/MultiplicityProbModelI.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISPartonModelPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"
#include "Utils/KineUtils.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchFx.h"

using std::ostringstream;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec() :
XSecAlgorithmI("genie::DISPartonModelPXSec")
{
  fInInitPhase = true;
}
//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec(string config) :
XSecAlgorithmI("genie::DISPartonModelPXSec", config)
{
  fInInitPhase = true;
}
//____________________________________________________________________________
DISPartonModelPXSec::~DISPartonModelPXSec()
{

}
//____________________________________________________________________________
double DISPartonModelPXSec::XSec(
                const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- Get kinematical & init-state parameters
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E     = init_state.GetProbeE(kRfStruckNucAtRest);
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double x     = kinematics.x();
  double y     = kinematics.y();

  LOG("DISPXSec", pDEBUG)  
   << "Computing d2xsec/dxdy @ E = " << E << ", x = " << x << ", y = " << y;

  //----- One of the xsec terms changes sign for antineutrinos
  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  //----- Calculate the DIS structure functions
  fDISSF.Calculate(interaction); 

  LOG("DISPXSec", pDEBUG) << fDISSF;

  //-- calculate auxiliary parameters
  double ml2  = ml    * ml;
  double ml4  = ml2   * ml2;
  double E2   = E     * E;
  double Gfac = (kGF*kGF*Mnuc*E) / kPi;

  //----- Build all dxsec/dxdy terms
  double term1 = y * ( x*y + ml2/(2*E*Mnuc) );
  double term2 = 1 - y - Mnuc*x*y/(2*E) - ml2/(4*E2);
  double term3 = sign * (x*y*(1-y/2) - y*ml2/(4*Mnuc*E));
  double term4 = x*y*ml2/(2*Mnuc*E) + ml4/(4*Mnuc2*E2);
  double term5 = -1.*ml2/(2*Mnuc*E);

  LOG("DISPXSec", pDEBUG)  
    << "\nd2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2 << ")*F2+(" 
                  << term3 << ")*F3+(" << term4 << ")*F4+(" << term5 << ")*F5";

  //----- Compute the differential cross section
  term1 *= fDISSF.F1();
  term2 *= fDISSF.F2();
  term3 *= fDISSF.F3();
  term4 *= fDISSF.F4();
  term5 *= fDISSF.F5();

  double xsec = Gfac*(term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec,0.);

  LOG("DISPXSec", pINFO)
        << "d2xsec/dxdy[FreeN] (E= " << E 
                      << ", x= " << x << ", y= " << y << ") = " << xsec;

  //----- If the DIS/RES joining scheme is enabled, modify the xsec accordingly
  if(fUsingDisResJoin) {
     double R = this->DISRESJoinSuppressionFactor(interaction);
     xsec*=R;
     LOG("DISPXSec", pINFO) << "d2xsec/dxdy[FreeN, D/R Join] " << xsec;
  }

  //----- The algorithm computes d^2xsec/dxdy
  //      Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  //----- If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- Compute nuclear cross section (simple scaling here, corrections must
  //      have been included in the structure functions)

  const Target & target = init_state.GetTarget();
  int nucpdgc = target.StruckNucleonPDGCode();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 
  xsec *= NNucl; 

  return xsec;
}
//____________________________________________________________________________
bool DISPartonModelPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool DISPartonModelPXSec::ValidKinematics(
                                       const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E     = init_state.GetProbeE(kRfStruckNucAtRest);
  double Mnuc  = init_state.GetTarget().StruckNucleonMass();
  double Mnuc2 = TMath::Power(Mnuc, 2);
  double x     = kinematics.x();
  double y     = kinematics.y();
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Sqrt(W2);
  double Q2    = 2*Mnuc*E*x*y;

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed
  Range1D_t rW  = utils::kinematics::KineRange (interaction, kKVW);
  Range1D_t rQ2 = utils::kinematics::KineRange (interaction, kKVQ2);

  bool in_range = utils::math::IsWithinLimits(Q2, rQ2)
                                      && utils::math::IsWithinLimits(W, rW);
  if(!in_range) {
       LOG("DISXSec", pDEBUG)
             << "*** (W = " << W
                           << ", Q2 = " << Q2 << " isn't in physical range";
       LOG("DISXSec", pDEBUG)
                << "Phys W: " << "[" << rW.min << ", " << rW.max << "] GeV";
       LOG("DISXSec", pDEBUG)
           << "Phys Q2: " << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
       return false;
  }
  return true;
}
//____________________________________________________________________________
double DISPartonModelPXSec::DISRESJoinSuppressionFactor(
                                                const Interaction * in) const
{
// Computes suppression factors for the DIS xsec under the used DIS/RES join
// scheme. Since this is a 'low-level' algorithm that is being called many
// times per generated event or computed cross section spline, the suppression
// factors would be cached to avoid calling the hadronic multiplicity model
// too often.
//
  double R=0, Ro=0;

  const double Wmin = kNeutronMass + kPionMass + 1E-3;

  const InitialState & ist = in->GetInitialState();
  const ProcessInfo &  pi  = in->GetProcessInfo();

  double E    = ist.GetProbeE(kRfStruckNucAtRest);
  double Mnuc = ist.GetTarget().StruckNucleonMass();
  double x    = in->GetKinematics().x(); 
  double y    = in->GetKinematics().y();
  double Wo   = utils::kinematics::XYtoW(E,Mnuc,x,y);

  if(!fUseCache) {
    // ** Compute the reduction factor at each call - no caching 
    //
    const TH1D & mprob = fMultProbModel->ProbabilityDistribution(in);
    R = mprob.Integral("width");
  }
  else {
    // ** Precompute/cache the reduction factors and then use the 
    // ** cache to evaluate these factors
    //

    //-- Access the cache branch. The branch key is formed as:
    //   algid/DIS-RES-Join/nu-pdg:N;hit-nuc-pdg:N/inttype

    Cache * cache = Cache::Instance();
    string algkey = this->Id().Key() + "/DIS-RES-Join";

    ostringstream ikey;
    ikey << "nu-pdgc:" << ist.GetProbePDGCode() 
         << ";hit-nuc-pdg:"<< ist.GetTarget().StruckNucleonPDGCode() << "/"
         << pi.InteractionTypeAsString();

    string key = cache->CacheBranchKey(algkey, ikey.str());

    CacheBranchFx * cbr =
          dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));

    //-- If it does not exist create a new one and cache DIS xsec suppression
    //   factors
    bool non_zero=false;
    if(!cbr) {
      LOG("DISXSec", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;

      cbr = new CacheBranchFx("DIS Suppr. Factors in DIS/RES Join Scheme");
      Interaction interaction(*in);

      const int kN   = 300;
      double WminSpl = Wmin;
      double WmaxSpl = fWcut + 0.1; // well into the area where scaling factor = 1
      double dW      = (WmaxSpl-WminSpl)/(kN-1);

      for(int i=0; i<kN; i++) {
        double W = WminSpl+i*dW;
        interaction.GetKinematicsPtr()->SetW(W);
        const TH1D & mprob = 
                    fMultProbModel->ProbabilityDistribution(&interaction);
        R = mprob.Integral("width");

        // make sure that it takes enough samples where it is non-zero:
        // modify the step and the sample counter once I've hit the first
        // non-zero value
        if(!non_zero && R>0) {
	  non_zero=true;
          WminSpl=W;
          i = 0;
          dW = (WmaxSpl-WminSpl)/(kN-1);
        }
        LOG("DISXSec", pNOTICE) 
	    << "Cached DIS XSec Suppr. factor (@ W=" << W << ") = " << R;

        cbr->AddValues(W,R);
      }
      cbr->CreateSpline();

      cache->AddCacheBranch(key, cbr);
      assert(cbr);
    } // cache data

    //-- get the reduction factor from the cache branch

    const CacheBranchFx & cache_branch = (*cbr);
    R = cache_branch(Wo);
  }

  //-- Now return the suppression factor

  if      (Wo > Wmin && Wo < fWcut-1E-2) Ro = R;
  else if (Wo <= Wmin)                   Ro = 0.0;
  else                                   Ro = 1.0;

  LOG("DISXSec", pDEBUG) 
            << "DIS/RES Join: DIS xsec suppr. (W=" << Wo << ") = " << Ro;
  return Ro;
}
//____________________________________________________________________________
void DISPartonModelPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPartonModelPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISPartonModelPXSec::LoadConfig(void)
{
  // Access global defaults to use in case of missing parameters
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fDISSFModel = 0;
  fDISSFModel = dynamic_cast<const DISStructureFuncModelI *> (
                                this->SubAlg("sf-alg-name", "sf-param-set"));
  assert(fDISSFModel);

  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm

  fUsingDisResJoin = fConfig->GetBoolDef(
		"use-dis-res-joining-scheme", gc->GetBool("UseDRJoinScheme"));
  fMultProbModel = 0;
  fWcut=0;

  if(fUsingDisResJoin) {
     fMultProbModel = 
         dynamic_cast<const MultiplicityProbModelI *> (this->SubAlg(
                 "multiplicity-prob-alg-name","multiplicity-prob-param-set"));
     assert(fMultProbModel);

     // Load Wcut determining the phase space area where the multiplicity prob.
     // scaling factors would be applied -if requested-
     fWcut = fConfig->GetDoubleDef("Wcut", gc->GetDouble("Wcut"));
  }

  //-- Caching the reduction factors used in the DIS/RES joing scheme?
  //   In normal event generation (1 config -> many calls) it is worth caching
  //   these suppression factors.
  //   Depending on the way this algorithm is used during event reweighting,
  //   precomputing (for all W's) & caching these factors might not be efficient.
  //   Here we provide the option to turn the caching off (default: on)

  fUseCache = fConfig->GetBoolDef("use-cache", true);

  //-- Since this method would be called every time the current algorithm is 
  //   reconfigured at run-time, remove all the data cached by this algorithm
  //   since they depend on the previous configuration

  if(!fInInitPhase) {
     Cache * cache = Cache::Instance();
     string keysubstr = this->Id().Key() + "/DIS-RES-Join";
     cache->RmMatchedCacheBranches(keysubstr);
  }
  fInInitPhase = false;
}
//____________________________________________________________________________

