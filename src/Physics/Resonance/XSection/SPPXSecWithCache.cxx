//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
  Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
  based on code of Costas Andreopoulos <c.andreopoulos \at cern.ch>
  University of Liverpool

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <algorithm>
#include <cassert>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/AdaptiveIntegratorMultiDim.h>

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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/Resonance/XSection/SPPXSecWithCache.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/XSecSplineList.h"

using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;

//____________________________________________________________________________
SPPXSecWithCache::SPPXSecWithCache() :
  XSecIntegratorI()
{

}
//____________________________________________________________________________
SPPXSecWithCache::SPPXSecWithCache(string nm) :
  XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
SPPXSecWithCache::SPPXSecWithCache(string nm,string conf):
  XSecIntegratorI(nm,conf)
{

}
//____________________________________________________________________________
SPPXSecWithCache::~SPPXSecWithCache()
{

}
//____________________________________________________________________________
void SPPXSecWithCache::CacheResExcitationXSec(const Interaction * in) const
{
  // Cache resonance neutrino production data from free nucleons

  Cache * cache = Cache::Instance();

  assert(fSinglePionProductionXSecModel);
  
  SppChannel_t spp_channel = SppChannel::FromInteraction(in);

  // Splines parameters are taken from Splines configuration
  XSecSplineList * xsl = XSecSplineList::Instance();

  // Compute the number of spline knots - use at least 10 knots per decade
  // && at least 40 knots in the full energy range
  const double Emin = xsl -> Emin() ;
  const double Emax = std::min( fEMax, xsl -> Emax() ) ; // here we ignore the run configuration 
                                                         // since we know that the splines have a plateau
  const int    nknots     = xsl -> NKnots() ;

  vector<double> E( nknots, 0. ) ; 
  
  int nu_code  = in->InitState().ProbePdg();
  int nuc_code = in->InitState().Tgt().HitNucPdg();
  int tgt_code = (nuc_code==kPdgProton) ? kPdgTgtFreeP : kPdgTgtFreeN;

  Interaction local_interaction(*in);
  local_interaction.InitStatePtr()->SetPdgs(tgt_code, nu_code);
  local_interaction.InitStatePtr()->TgtPtr()->SetHitNucPdg(nuc_code);

  InteractionType_t wkcur = local_interaction.ProcInfo().InteractionTypeId();

  // Get a unique cache branch name
  string key = this->CacheBranchName(spp_channel, wkcur, nu_code);

  // Make sure the cache branch does not already exists
  CacheBranchFx * cache_branch =
    dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
  assert(!cache_branch);
  
  // Create the new cache branch
  LOG("SPPCache", pNOTICE)
    << "\n ** Creating cache branch - key = " << key;
  cache_branch = new CacheBranchFx("ResSPP XSec");
  cache->AddCacheBranch(key, cache_branch);
  assert(cache_branch);
  
  const KPhaseSpace& kps = in->PhaseSpace();
    
  double Ethr = kps.Threshold_SPP_iso();
  LOG("SPPCache", pNOTICE) << "E threshold = " << Ethr;

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
  double dEa = (TMath::Log10(Emax) - TMath::Log10(E0)) /(nka-1);
  for(int i=0; i<nka; i++) {
    E[i+nkb] = TMath::Power(10., TMath::Log10(E0) + i*dEa);
  }
  
  // Compute cross sections at the given set of energies
  for(int ie=0; ie<nknots; ie++) {
    double xsec = 0.;
    double Ev   = E[ie];
    
    TLorentzVector p4(0., 0., Ev, Ev);
    local_interaction.InitStatePtr()->SetProbeP4(p4);
    
    if(Ev>Ethr+kASmallNum) {
      
      LOG("SPPCache", pINFO)
      << "*** Integrating d^3 XSec/dWdQ^2dCosTheta for Ch: "
      << SppChannel::AsString(spp_channel) << " at Ev = " << Ev;
      
      utils::gsl::d3XSecMK_dWQ2CosTheta_E func(fSinglePionProductionXSecModel, & local_interaction, fWcut ) ; 
      ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
      ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
      if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) 
      {
        ROOT::Math::AdaptiveIntegratorMultiDim * cast = dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
        assert(cast);
        cast->SetMinPts(fGSLMinEval);
      }
      ig.SetFunction(func);
      double kine_min[3] = { 0., 0., 0.};
      double kine_max[3] = { 1., 1., 1.};
      xsec = ig.Integral(kine_min, kine_max)*(1E-38 * units::cm2);;
    } 
    else 
      LOG("SPPCache", pINFO) << "** Below threshold E = " << Ev << " <= " << Ethr;
    
    cache_branch->AddValues(Ev,xsec);
    
    
    string nc_nuc   = ((pdg::IsNeutrino(nu_code)) ? ";v:" : ";vb:");
    
    
    SLOG("SPPCache", pNOTICE)
      << "SPP XSec (Ch:" << SppChannel::AsString(spp_channel) << nc_nuc  << nu_code
      << ", E="<< Ev << ") = "<< xsec << " x 1E-38 cm^2";
  }//spline knots
  
  // Build the spline
  cache_branch->CreateSpline();
    
}
//____________________________________________________________________________
string SPPXSecWithCache::CacheBranchName(
                    SppChannel_t spp_channel, InteractionType_t it, int nupdgc) const
{
  // Build a unique name for the cache branch

  Cache * cache = Cache::Instance();
  string spp_channel_name = SppChannel::AsString(spp_channel);
  string it_name  = InteractionType::AsString(it);
  string nc_nuc   = ((pdg::IsNeutrino(nupdgc)) ? ";v:" : ";vb:");
  
  ostringstream intk;
  intk << "ResSPPXSec/Ch:" << spp_channel_name << nc_nuc  << nupdgc
       << ";int:" << it_name;

  string algkey = fSinglePionProductionXSecModel->Id().Key();
  string ikey   = intk.str();
  string key    = cache->CacheBranchKey(algkey, ikey);

  return key;
}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::d3XSecMK_dWQ2CosTheta_E(
                                const XSecAlgorithmI * m, const Interaction * interaction, double  wcut) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fWcut(wcut)
{

  isZero = false;
  fInteraction = const_cast<Interaction*>(interaction);
  fInteraction->SetBit(kISkipProcessChk);
  fInteraction->SetBit(kISkipKinematicChk);
  
  kps = fInteraction->PhaseSpacePtr();
  
  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  double Enu = init_state.ProbeE(kRfHitNucRest);


  if (Enu < kps->Threshold_SPP_iso())
  {
    isZero = true;
    return;
  }
  
  Wl  = kps->WLim_SPP_iso();
  if (fWcut >= Wl.min)
    Wl.max = TMath::Min(fWcut,Wl.max);
  
  
}
genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::~d3XSecMK_dWQ2CosTheta_E()
{

}
unsigned int genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::DoEval(const double * xin) const
{
  // outputs:
  //   differential cross section [1/GeV^3] for Resonance single pion production production
  //

  if (isZero) return 0.;
  
  double W2  = Wl.min*Wl.min + (Wl.max*Wl.max - Wl.min*Wl.min)*xin[0];
  fInteraction->KinePtr()->SetW(TMath::Sqrt(W2));
   
  Range1D_t Q2l = kps->Q2Lim_W_SPP_iso(); 
   
  double sqrt_Q2 = TMath::Sqrt(Q2l.min) + ( TMath::Sqrt(Q2l.max) - TMath::Sqrt(Q2l.min) )*xin[1];
  fInteraction->KinePtr()->SetQ2(sqrt_Q2*sqrt_Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[2]); //CosTheta
  
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpfE)*sqrt_Q2*(Wl.max*Wl.max - Wl.min*Wl.min)*(TMath::Sqrt(Q2l.max) - TMath::Sqrt(Q2l.min))*2/TMath::Sqrt(W2);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E(fModel,fInteraction,fWcut);
}
