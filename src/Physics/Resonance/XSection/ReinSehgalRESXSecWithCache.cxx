//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - March 09, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 18, 2008 - CA
   Simplify the way free nucleon channels (for which we cache cross sections)
   are built from the input interaction
 @ Jan 19, 2008 - CA
   Modify the way knots are distributed in the cached free nucleon resonance
   neutrino production splines so that the energy threshold is treated more
   accurately (see also XSecSplineList.cxx).
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Physics/Resonance/XSection/ReinSehgalRESXSecWithCache.h"

using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//using namespace genie::units;

//____________________________________________________________________________
ReinSehgalRESXSecWithCache::ReinSehgalRESXSecWithCache() :
XSecIntegratorI()
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCache::ReinSehgalRESXSecWithCache(string nm) :
XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCache::ReinSehgalRESXSecWithCache(string nm,string conf):
XSecIntegratorI(nm,conf)
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCache::~ReinSehgalRESXSecWithCache()
{

}
//____________________________________________________________________________
void ReinSehgalRESXSecWithCache::CacheResExcitationXSec(
                                                 const Interaction * in) const
{
// Cache resonance neutrino production data from free nucleons

  Cache * cache = Cache::Instance();

  assert(fSingleResXSecModel);
//  assert(fIntegrator);

  // Compute the number of spline knots - use at least 10 knots per decade 
  // && at least 40 knots in the full energy range
  const double Emin       = 0.01;
  const int    nknots_min = (int) (10*(TMath::Log(fEMax)-TMath::Log(Emin))); 
  const int    nknots     = TMath::Max(40, nknots_min); 
  double * E = new double[nknots]; // knot 'x'

  TLorentzVector p4(0,0,0,0);

  int nu_code  = in->InitState().ProbePdg();
  int nuc_code = in->InitState().Tgt().HitNucPdg();
  int tgt_code = (nuc_code==kPdgProton) ? kPdgTgtFreeP : kPdgTgtFreeN;

  Interaction * interaction = new Interaction(*in);
  interaction->InitStatePtr()->SetPdgs(tgt_code, nu_code);
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(nuc_code);

  InteractionType_t wkcur = interaction->ProcInfo().InteractionTypeId();
  unsigned int nres = fResList.NResonances();
  for(unsigned int ires = 0; ires < nres; ires++) {

         // Get next resonance from the resonance list
         Resonance_t res = fResList.ResonanceId(ires);

         interaction->ExclTagPtr()->SetResonance(res);

         // Get a unique cache branch name
         string key = this->CacheBranchName(res, wkcur, nu_code, nuc_code);

         // Make sure the cache branch does not already exists
         CacheBranchFx * cache_branch =
             dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
         assert(!cache_branch);

         // Create the new cache branch
         LOG("ReinSehgalResC", pNOTICE) 
                        << "\n ** Creating cache branch - key = " << key;
         cache_branch = new CacheBranchFx("RES Excitation XSec");
         cache->AddCacheBranch(key, cache_branch);
         assert(cache_branch);

         const KPhaseSpace & kps = interaction->PhaseSpace();
         double Ethr = kps.Threshold();
         LOG("ReinSehgalResC", pNOTICE) << "E threshold = " << Ethr;

         // Distribute the knots in the energy range as is being done in the
         // XSecSplineList so that the energy threshold is treated correctly
         // in the spline - see comments there in.
         int nkb = (Ethr>Emin) ? 5 : 0; // number of knots <  threshold
         int nka = nknots-nkb;          // number of knots >= threshold
         // knots < energy threshold
         double dEb =  (Ethr>Emin) ? (Ethr - Emin) / nkb : 0;
         for(int i=0; i<nkb; i++) {
            E[i] = Emin + i*dEb;
         }
         // knots >= energy threshold
         double E0  = TMath::Max(Ethr,Emin);
         double dEa = (TMath::Log10(fEMax) - TMath::Log10(E0)) /(nka-1);
         for(int i=0; i<nka; i++) {
            E[i+nkb] = TMath::Power(10., TMath::Log10(E0) + i * dEa);
         }

         // Compute cross sections at the given set of energies
         for(int ie=0; ie<nknots; ie++) {
             double xsec = 0.;
             double Ev   = E[ie];
             p4.SetPxPyPzE(0,0,Ev,Ev);
             interaction->InitStatePtr()->SetProbeP4(p4);

             if(Ev>Ethr+kASmallNum) {
               // Get W integration range and the wider possible Q2 range 
               // (for all W)
               Range1D_t rW  = kps.Limits(kKVW);
               Range1D_t rQ2 = kps.Limits(kKVQ2);

               LOG("ReinSehgalResC", pINFO) 
	         << "*** Integrating d^2 XSec/dWdQ^2 for R: " 
     	                 << utils::res::AsString(res) << " at Ev = " << Ev;
               LOG("ReinSehgalResC", pINFO) 
                                << "{W}   = " << rW.min  << ", " << rW.max;
      	       LOG("ReinSehgalResC", pINFO) 
                               << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

               if(rW.max<rW.min || rQ2.max<rQ2.min || rW.min<0 || rQ2.min<0) {
     	          LOG("ReinSehgalResC", pINFO) 
                              << "** Not allowed kinematically, xsec=0";
               } else {

                  ROOT::Math::IBaseFunctionMultiDim * func = 
                      new utils::gsl::d2XSec_dWdQ2_E(fSingleResXSecModel, interaction);
                  ROOT::Math::IntegrationMultiDim::Type ig_type = 
                      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
                  ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);   
                  ig.SetFunction(*func);
                  double kine_min[2] = { rW.min, rQ2.min };
                  double kine_max[2] = { rW.max, rQ2.max };
                  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
                  delete func;
               }
             } else {
                 LOG("ReinSehgalResC", pINFO) 
 		       << "** Below threshold E = " << Ev << " <= " << Ethr;
             }
             cache_branch->AddValues(Ev,xsec);
             SLOG("ReinSehgalResC", pNOTICE) 
               << "RES XSec (R:" << utils::res::AsString(res)
    	       << ", E="<< Ev << ") = "<< xsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";
         }//spline knots

         // Build the spline
         cache_branch->CreateSpline();
  }//ires

  delete [] E;
  delete interaction;
}
//____________________________________________________________________________
string ReinSehgalRESXSecWithCache::CacheBranchName(
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
