//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 09, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "ReinSeghal/ReinSeghalRESXSecWithCache.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchFx.h"

using std::ostringstream;
using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache() :
XSecAlgorithmI()
{

}
//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache(string name) :
XSecAlgorithmI(name)
{

}
//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache(
                                                 string name, string config) :
XSecAlgorithmI(name,config)
{

}
//____________________________________________________________________________
ReinSeghalRESXSecWithCache::~ReinSeghalRESXSecWithCache()
{

}
//____________________________________________________________________________
void ReinSeghalRESXSecWithCache::CacheResExcitationXSec(
                                                 const Interaction * in) const
{
// Cache resonance neutrino production data from free nucleons

  Cache * cache = Cache::Instance();

  assert(fSingleResXSecModel);
  assert(fIntegrator);

  Interaction * interaction = new Interaction(*in);
  interaction->TestBit(kIAssumeFreeNucleon);

  const int kNSt = 3;
  const InteractionType_t it[kNSt] = {kIntWeakCC, kIntWeakNC, kIntWeakNC };
  const int iN[kNSt] = {kPdgProton, kPdgProton, kPdgNeutron};

  int nu_code  = interaction->GetInitialState().GetProbePDGCode();

  // at the splines use at least 10 knots per decade but at least 40 knots
  // in the full energy range
  const double kEmin         = 0.120; // xsec_res(Emin) = 0
  const double kLogEmin      = TMath::Log(kEmin);
  const double kLogEmax      = TMath::Log(fEMax);
  const int    kMinNKnots    = (int) (10*(kLogEmax-kLogEmin)); 
  const int    kNSplineKnots = TMath::Max(40, kMinNKnots); 
  const double kdLogE        = (kLogEmax-kLogEmin)/(kNSplineKnots-1);

  unsigned int nres = fResList.NResonances();
  for(int i=0; i<kNSt; i++) {

    interaction->GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(iN[i]);
    interaction->GetProcessInfoPtr()->Set(kScResonant, it[i]);

    for(unsigned int ires = 0; ires < nres; ires++) {

      //-- Get next resonance from the resonance list
      Resonance_t res = fResList.ResonanceId(ires);

      interaction->GetExclusiveTagPtr()->SetResonance(res);

      //-- Get a unique cache branch name
      string key = this->CacheBranchName(res, it[i], nu_code, iN[i]);

      //-- Make sure the cache branch does not already exists
      CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
      assert(!cache_branch);

      //-- Create the new cache branch
      LOG("ReinSeghalResC", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;
      cache_branch = new CacheBranchFx("RES Excitation XSec");
      cache->AddCacheBranch(key, cache_branch);
      assert(cache_branch);

      TLorentzVector p4(0,0,0,0);

      for(int ie=0; ie<kNSplineKnots; ie++) {

        double Ev = TMath::Exp(kLogEmin + ie*kdLogE);
        p4.SetPxPyPzE(0,0,Ev,Ev);
        interaction->GetInitialStatePtr()->SetProbeP4(p4);

        // Get W integration range
        Range1D_t rW = this->WRange(interaction);
        // Get the wider possible Q2 range for the input W range
        Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction, rW);

        LOG("ReinSeghalResC", pINFO) 
	  << "*** Integrating d^2 XSec/dWdQ^2 for R: " 
     	      << utils::res::AsString(res) << " at Ev = " << Ev;
        LOG("ReinSeghalResC", pINFO) 
                     << "{W}   = " << rW.min  << ", " << rW.max;
	LOG("ReinSeghalResC", pINFO) 
                    << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

        double xsec = 0;

        if(rW.max<rW.min || rQ2.max<rQ2.min || rW.min<0 || rQ2.min<0) {
	  LOG("ReinSeghalResC", pINFO) 
                              << "** Not allowed kinematically, xsec=0";
        } else {
          GXSecFunc * func = new Integrand_D2XSec_DWDQ2_E(
                                      fSingleResXSecModel, interaction);
          func->SetParam(0,"W",  rW);
          func->SetParam(1,"Q2", rQ2);
          xsec = fIntegrator->Integrate(*func);
          delete func;
        }
        cache_branch->AddValues(Ev,xsec);
        LOG("ReinSeghalResC", pNOTICE) 
            << "RES XSec (R:" << utils::res::AsString(res)
                                       << ", E=" << Ev << ")= " << xsec;
      }//ie

      cache_branch->CreateSpline();

    }//ires
  }//i

  delete interaction;
}
//____________________________________________________________________________
string ReinSeghalRESXSecWithCache::CacheBranchName(
     Resonance_t res, InteractionType_t it, int nupdgc, int nucleonpdgc) const
{
// Build a unique name for the cache branch

  Cache * cache = Cache::Instance();
  string res_name = utils::res::AsString(res);
  string it_name  = InteractionType::AsString(it);
  string nc_nuc   = "";
  if(it == kIntWeakNC) { nc_nuc = ((nucleonpdgc==kPdgProton) ? "p" : "n"); }

  ostringstream intk;
  intk << "ResExcitationXSec/R:" << res_name << ";nu:"  << nupdgc
           << ";int:" << it_name << nc_nuc;
      
  string algkey = fSingleResXSecModel->Id().Key();
  string ikey   = intk.str();
  string key    = cache->CacheBranchKey(algkey, ikey);

  return key;
}
//____________________________________________________________________________
Range1D_t ReinSeghalRESXSecWithCache::WRange(
                                      const Interaction * interaction) const
{
  //-- Get the physically allowed W range for this interaction and allow the
  //   user inputs (if any) to narrow it
  Range1D_t rW = utils::kinematics::KineRange(interaction, kKVW);
  LOG("ReinSeghalResC", pDEBUG)
      << "Physical W range: " << "[" << rW.min << ", " << rW.max << "] GeV";
  // apply user cuts
  utils::kinematics::ApplyCutsToKineLimits(rW, fWminCut,  fWmaxCut );
  LOG("ReinSeghalResC", pDEBUG)
       << "Physical & User W range: "
                              << "[" << rW.min << ", " << rW.max << "] GeV";
  return rW;
}
//___________________________________________________________________________
Range1D_t ReinSeghalRESXSecWithCache::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed Q2 range for this interaction and allow the
  //   user inputs (if any) to narrow it
  Range1D_t rQ2 = utils::kinematics::KineRange(interaction, kKVQ2); 
  LOG("ReinSeghalResC", pDEBUG) << "Physical Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  // apply user cuts
  utils::kinematics::ApplyCutsToKineLimits(rQ2, fQ2minCut, fQ2maxCut);
  LOG("ReinSeghalResC", pDEBUG)
       << "Physical & User Q2 range: "
                         << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
  return rQ2;
}
//___________________________________________________________________________


