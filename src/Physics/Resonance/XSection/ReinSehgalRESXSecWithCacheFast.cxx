//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin <kakorin@jinr.ru>
 Joint Institute for Nuclear Research

 based on code of
 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <sstream>
#include <cassert>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/Resonance/XSection/ReinSehgalRESXSecWithCacheFast.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Utils/Range1.h"







using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
//using namespace genie::units;

//____________________________________________________________________________
ReinSehgalRESXSecWithCacheFast::ReinSehgalRESXSecWithCacheFast() :
XSecIntegratorI()
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCacheFast::ReinSehgalRESXSecWithCacheFast(string nm) :
XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCacheFast::ReinSehgalRESXSecWithCacheFast(string nm,string conf):
XSecIntegratorI(nm,conf)
{

}
//____________________________________________________________________________
ReinSehgalRESXSecWithCacheFast::~ReinSehgalRESXSecWithCacheFast()
{

}
//____________________________________________________________________________
void ReinSehgalRESXSecWithCacheFast::CacheResExcitationXSec(
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
  const int    nknots     = TMath::Max(100, nknots_min);
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
         LOG("ReinSehgalResCF", pNOTICE)
                        << "\n ** Creating cache branch - key = " << key;
         cache_branch = new CacheBranchFx("RES Excitation XSec");
         cache->AddCacheBranch(key, cache_branch);
         assert(cache_branch);

         const KPhaseSpace & kps = interaction->PhaseSpace();
         double Ethr = kps.Threshold();
         LOG("ReinSehgalResCF", pNOTICE) << "E threshold = " << Ethr;

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
             // Get integration ranges
                                Range1D_t rW  = Range1D_t(0.0,1.0);
                                Range1D_t rQ2 = Range1D_t(0.0,1.0);

               LOG("ReinSehgalResCF", pINFO)
                 << "*** Integrating d^2 XSec/dWdQ^2 for R: "
                         << utils::res::AsString(res) << " at Ev = " << Ev;
               LOG("ReinSehgalResCF", pINFO)
                                << "{W}   = " << rW.min  << ", " << rW.max;
               LOG("ReinSehgalResCF", pINFO)
                               << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

               if(rW.max<rW.min || rQ2.max<rQ2.min || rW.min<0 || rQ2.min<0) {
                  LOG("ReinSehgalResCF", pINFO)
                              << "** Not allowed kinematically, xsec=0";
               } else {
                                  ROOT::Math::IBaseFunctionMultiDim * func= new utils::gsl::d2XSecRESFast_dWQ2_E(fSingleResXSecModel, interaction);
                  ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
                  ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
                  ig.SetFunction(*func);
                  double kine_min[2] = { rW.min, rQ2.min };
                  double kine_max[2] = { rW.max, rQ2.max };
                  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
                  delete func;
               }
             } else {
                 LOG("ReinSehgalResCF", pINFO)
                       << "** Below threshold E = " << Ev << " <= " << Ethr;
             }
             cache_branch->AddValues(Ev,xsec);
             SLOG("ReinSehgalResCF", pNOTICE)
               << "RES XSec (R:" << utils::res::AsString(res)
               << ", E="<< Ev << ") = "<< xsec/(1E-38 *genie::units::cm2)
               << " x 1E-38 cm^2";
         }//spline knots

         // Build the spline
         cache_branch->CreateSpline();
  }//ires

  delete [] E;
  delete interaction;
}
//____________________________________________________________________________
string ReinSehgalRESXSecWithCacheFast::CacheBranchName(
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
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d2XSecRESFast_dWQ2_E::d2XSecRESFast_dWQ2_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
        kps = fInteraction->PhaseSpacePtr();
        Range1D_t Wl  = kps->WLim();
        fWmin = Wl.min;
        fWmax = Wl.max;
        Registry fConfig = (const_cast<XSecAlgorithmI *>(fModel))->GetConfig();
        bool fUsingDisResJoin = fConfig.GetBool("UseDRJoinScheme");
        double fWcut = 999999;
        if(fUsingDisResJoin)
        {
                fWcut = fConfig.GetDouble("Wcut");
        }
        fWmax=TMath::Min(fWcut, fWmax);
        if (fWcut<fWmin)
                isfWcutLessfWmin=true;
        else
                isfWcutLessfWmin=false;
        bool fNormBW = fConfig.GetBoolDef("BreitWignerNorm", true);
        if (fNormBW)
        {
                double fN2ResMaxNWidths = fConfig.GetDoubleDef("MaxNWidthForN2Res", 2.0);
                double fN0ResMaxNWidths = fConfig.GetDoubleDef("MaxNWidthForN0Res", 6.0);
                double fGNResMaxNWidths = fConfig.GetDoubleDef("MaxNWidthForGNRes", 4.0);
                Resonance_t resonance = fInteraction->ExclTag().Resonance();
                int    IR  = utils::res::ResonanceIndex    (resonance);
                double MR  = utils::res::Mass              (resonance);
                double WR  = utils::res::Width             (resonance);
                if (IR==0)
                        fWcut = MR + fN0ResMaxNWidths * WR;
                else if (IR==2)
                        fWcut = MR + fN2ResMaxNWidths * WR;
                else
                        fWcut = MR + fGNResMaxNWidths * WR;
                fWmax=TMath::Min(fWcut, fWmax);
                if (fWcut<fWmin)
                        isfWcutLessfWmin=true;
        }

}
genie::utils::gsl::d2XSecRESFast_dWQ2_E::~d2XSecRESFast_dWQ2_E()
{

}
unsigned int genie::utils::gsl::d2XSecRESFast_dWQ2_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::d2XSecRESFast_dWQ2_E::DoEval(const double * xin) const
{
// inputs:
//    W
//    Q2
// outputs:
//   differential cross section [10^-38 cm^2/GeV^3] for Resonance production
//

  if (isfWcutLessfWmin)
        return 0;
  double W  = fWmin+(fWmax-fWmin)*xin[0];
  fInteraction->KinePtr()->SetW(W);
  Range1D_t Q2l = kps->Q2Lim_W();
  if (Q2l.min<0 || Q2l.max<0)
        return 0.0;
  double Q2 = Q2l.min+(Q2l.max-Q2l.min)*xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE)*(fWmax-fWmin)*(Q2l.max-Q2l.min);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2XSecRESFast_dWQ2_E::Clone() const
{
  return
    new genie::utils::gsl::d2XSecRESFast_dWQ2_E(fModel,fInteraction);
}
//____________________________________________________________________________
