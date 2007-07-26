//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - March 09, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
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
using namespace genie::units;

//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache() :
XSecIntegratorI()
{

}
//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache(string nm) :
XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
ReinSeghalRESXSecWithCache::ReinSeghalRESXSecWithCache(string nm,string conf):
XSecIntegratorI(nm,conf)
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

  const int kNTgt = 2;
  const int kNWkC = 2;
  const int               tgtc[kNTgt] = { kPdgTgtFreeP, kPdgTgtFreeN };
  const InteractionType_t wkcc[kNWkC] = { kIntWeakCC, kIntWeakNC     };

  // at the splines use at least 10 knots per decade but at least 40 knots
  // in the full energy range
  const double kEmin         = 0.120; // xsec_res(Emin) = 0
  const double kLogEmin      = TMath::Log(kEmin);
  const double kLogEmax      = TMath::Log(fEMax);
  const int    kMinNKnots    = (int) (10*(kLogEmax-kLogEmin)); 
  const int    kNSplineKnots = TMath::Max(40, kMinNKnots); 
  const double kdLogE        = (kLogEmax-kLogEmin)/(kNSplineKnots-1);

  TLorentzVector p4(0,0,0,0);

  int nu_code = in->InitState().ProbePdg();

  Interaction * interaction = new Interaction(*in);
  interaction->TestBit(kIAssumeFreeNucleon);

  InitialState * init_state = interaction->InitStatePtr();
  ProcessInfo *  proc_info  = interaction->ProcInfoPtr();

  unsigned int nres = fResList.NResonances();

  for(int iwkc=0; iwkc<kNWkC; iwkc++) {
    for(int itgt=0; itgt<kNTgt; itgt++) {

      init_state -> SetPdgs(tgtc[itgt], nu_code);
      proc_info  -> Set(kScResonant, wkcc[iwkc]);

      for(unsigned int ires = 0; ires < nres; ires++) {

         //-- Get next resonance from the resonance list
         Resonance_t res = fResList.ResonanceId(ires);

         interaction->ExclTagPtr()->SetResonance(res);

         //-- Get a unique cache branch name
         int nuc_code = init_state->Tgt().HitNucPdg();
         string key = this->CacheBranchName(
                                   res, wkcc[iwkc], nu_code, nuc_code);

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

         const KPhaseSpace & kps = interaction->PhaseSpace();
         double EvThr = kps.Threshold();
         //double EvThr = interaction->EnergyThreshold();
         LOG("ReinSeghalResC", pNOTICE) << "E threshold = " << EvThr;

         for(int ie=0; ie<kNSplineKnots; ie++) {

             double xsec = 0.;
             double Ev   = TMath::Exp(kLogEmin + ie*kdLogE);
             p4.SetPxPyPzE(0,0,Ev,Ev);
             interaction->InitStatePtr()->SetProbeP4(p4);

             if(Ev>EvThr) {
               // Get W integration range and the wider possible Q2 range 
               // (for all W)
               Range1D_t rW  = kps.Limits(kKVW);
               Range1D_t rQ2 = kps.Limits(kKVQ2);

               LOG("ReinSeghalResC", pINFO) 
	         << "*** Integrating d^2 XSec/dWdQ^2 for R: " 
     	                 << utils::res::AsString(res) << " at Ev = " << Ev;
               LOG("ReinSeghalResC", pINFO) 
                                << "{W}   = " << rW.min  << ", " << rW.max;
      	       LOG("ReinSeghalResC", pINFO) 
                               << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

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
             } else {
                 LOG("ReinSeghalResC", pINFO) 
 		       << "** Below threshold E = " << Ev << " <= " << EvThr;
             }
             cache_branch->AddValues(Ev,xsec);
             SLOG("ReinSeghalResC", pNOTICE) 
               << "RES XSec (R:" << utils::res::AsString(res)
    	       << ", E="<< Ev << ") = "<< xsec/(1E-38 *cm2)<< " x 1E-38 cm^2";
         }//spline knots

         cache_branch->CreateSpline();

      }//ires
    }//hit nucleon
  }//weak current

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
  string nc_nuc   = ((nucleonpdgc==kPdgProton) ? "p" : "n"); 

  ostringstream intk;
  intk << "ResExcitationXSec/R:" << res_name << ";nu:"  << nupdgc
           << ";int:" << it_name << nc_nuc;
      
  string algkey = fSingleResXSecModel->Id().Key();
  string ikey   = intk.str();
  string key    = cache->CacheBranchKey(algkey, ikey);

  return key;
}
//____________________________________________________________________________
