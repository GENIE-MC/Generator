//____________________________________________________________________________
/*
  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
  Konstantin Kuzmin <kkuzmin@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
  Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
  based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
  University of Liverpool & STFC Rutherford Appleton Lab

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <sstream>
#include <algorithm>
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
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/Resonance/XSection/MKSPPXSecWithCache.h"
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
MKSPPXSecWithCache::MKSPPXSecWithCache() :
  XSecIntegratorI()
{

}
//____________________________________________________________________________
MKSPPXSecWithCache::MKSPPXSecWithCache(string nm) :
  XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
MKSPPXSecWithCache::MKSPPXSecWithCache(string nm,string conf):
  XSecIntegratorI(nm,conf)
{

}
//____________________________________________________________________________
MKSPPXSecWithCache::~MKSPPXSecWithCache()
{

}
//____________________________________________________________________________
void MKSPPXSecWithCache::CacheResExcitationXSec(const Interaction * in) const
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
  LOG("MKSPPCache", pNOTICE)
    << "\n ** Creating cache branch - key = " << key;
  cache_branch = new CacheBranchFx("ResSPP XSec");
  cache->AddCacheBranch(key, cache_branch);
  assert(cache_branch);
  
  const KPhaseSpace& kps = in->PhaseSpace();
    
  double Ethr = kps.Threshold_RSPP();
  LOG("MKSPPCache", pNOTICE) << "E threshold = " << Ethr;

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
      
      LOG("MKSPPCache", pINFO)
	<< "*** Integrating d^4 XSec/dWdQ^2dCosThetadPhi for Ch: "
	<< SppChannel::AsString(spp_channel) << " at Ev = " << Ev;
      
      utils::gsl::d3XSecMK_dWQ2CosTheta_E func(fSinglePionProductionXSecModel, & local_interaction, fWcut ) ; 
      ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
      ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
      ig.SetFunction(func);
      double kine_min[3] = { 0., 0., 0.};
      double kine_max[3] = { 1., 1., 1.};
      xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
      
    } 
    else 
      LOG("MKSPPCache", pINFO) << "** Below threshold E = " << Ev << " <= " << Ethr;
    
    cache_branch->AddValues(Ev,xsec);
    
    
    string nc_nuc   = ((pdg::IsNeutrino(nu_code)) ? ";v:" : ";vb:");
    
    
    SLOG("MKSPPCache", pNOTICE)
      << "ResSPP XSec (Ch:" << SppChannel::AsString(spp_channel) << nc_nuc  << nu_code
      << ", E="<< Ev << ") = "<< xsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";
  }//spline knots
  
  // Build the spline
  cache_branch->CreateSpline();
    
}
//____________________________________________________________________________
string MKSPPXSecWithCache::CacheBranchName(
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


  if (Enu < kps->Threshold_RSPP())
  {
    isZero = true;
    return;
  }
  
  Wl  = kps->WLim_RSPP();
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
  //   differential cross section [10^-38 cm^2/GeV^3] for Resonance single pion production production
  //

  if (isZero) return 0.;
  
  double W  = Wl.min + (Wl.max - Wl.min)*xin[0];
  fInteraction->KinePtr()->SetW(W);
   
  Range1D_t Q2l = kps->Q2Lim_W_RSPP(); 
   
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[2]); //CosTheta
  
  
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpfE)*(Wl.max-Wl.min)*(Q2l.max-Q2l.min)*2;
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSecMK_dWQ2CosTheta_E(fModel,fInteraction,fWcut);
}
