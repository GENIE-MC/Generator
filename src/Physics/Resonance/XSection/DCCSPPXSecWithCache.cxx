//____________________________________________________________________________
/*
  Copyright (c) 2003-2022, GENIE Neutrino MC Generator Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
  based on code of Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
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
#include "Physics/Resonance/XSection/DCCSPPXSecWithCache.h"
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
DCCSPPXSecWithCache::DCCSPPXSecWithCache() :
  XSecIntegratorI()
{

}
//____________________________________________________________________________
DCCSPPXSecWithCache::DCCSPPXSecWithCache(string nm) :
  XSecIntegratorI(nm)
{

}
//____________________________________________________________________________
DCCSPPXSecWithCache::DCCSPPXSecWithCache(string nm,string conf):
  XSecIntegratorI(nm,conf)
{

}
//____________________________________________________________________________
DCCSPPXSecWithCache::~DCCSPPXSecWithCache()
{

}
//____________________________________________________________________________
void DCCSPPXSecWithCache::CacheResExcitationXSec(const Interaction * in) const
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
  
  int nu_code         = in->InitState().ProbePdg();
  int nuc_code        = in->InitState().Tgt().HitNucPdg();
  int tgt_code        = (nuc_code==kPdgProton) ? kPdgTgtFreeP : kPdgTgtFreeN;
  int probe_helicity  = in->InitState().ProbeHelicity();

  Interaction local_interaction(*in);
  local_interaction.InitStatePtr()->SetPdgs(tgt_code, nu_code);
  local_interaction.InitStatePtr()->TgtPtr()->SetHitNucPdg(nuc_code);

  InteractionType_t wkcur = local_interaction.ProcInfo().InteractionTypeId();

  // Get a unique cache branch name
  string key = this->CacheBranchName(spp_channel, wkcur, nu_code, probe_helicity);

  // Make sure the cache branch does not already exists
  CacheBranchFx * cache_branch =
    dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
  assert(!cache_branch);
  
  // Create the new cache branch
  LOG("DCCSPPCache", pNOTICE)
    << "\n ** Creating cache branch - key = " << key;
  cache_branch = new CacheBranchFx("ResSPP XSec");
  cache->AddCacheBranch(key, cache_branch);
  assert(cache_branch);
  
  const KPhaseSpace& kps = in->PhaseSpace();
    
  double Ethr = kps.Threshold();
  LOG("DCCSPPCache", pNOTICE) << "E threshold = " << Ethr;

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
      
      LOG("DCCSPPCache", pINFO)
    << "*** Integrating d^3 XSec/dWdQ^2dCosTheta for Ch: "
    << SppChannel::AsString(spp_channel) << " at Ev = " << Ev;
      
      utils::gsl::d3XSecSPP_dWQ2CosTheta_E func(fSinglePionProductionXSecModel, & local_interaction, fWcut ) ; 
      ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
      ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
      ig.SetFunction(func);
      double kine_min[3] = { 0., 0., 0.};
      double kine_max[3] = { 1., 1., 1.};
      xsec = ig.Integral(kine_min, kine_max);
      
    } 
    else 
      LOG("DCCSPPCache", pINFO) << "** Below threshold E = " << Ev << " <= " << Ethr;
    
    cache_branch->AddValues(Ev,xsec);
    
    string nc_nuc   = this->ProbeAsString(nu_code, probe_helicity);
    
    
    SLOG("DCCSPPCache", pNOTICE)
      << "ResSPP XSec (Ch:" << SppChannel::AsString(spp_channel) << nc_nuc
      << ", E="<< Ev << ") = "<< xsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";
  }//spline knots
  
  // Build the spline
  cache_branch->CreateSpline();
    
}
//____________________________________________________________________________
string DCCSPPXSecWithCache::CacheBranchName(
                                  SppChannel_t spp_channel, InteractionType_t it, int nupdgc, int helicity) const
{
  // Build a unique name for the cache branch

  Cache * cache = Cache::Instance();
  string spp_channel_name = SppChannel::AsString(spp_channel);
  string it_name  = InteractionType::AsString(it);
  string nc_nuc = this->ProbeAsString(nupdgc, helicity);
  
  ostringstream intk;
  intk << "ResSPPXSec/Ch:" << spp_channel_name << nc_nuc << ";int:" << it_name;

  string algkey = fSinglePionProductionXSecModel->Id().Key();
  string ikey   = intk.str();
  string key    = cache->CacheBranchKey(algkey, ikey);

  return key;
}
//____________________________________________________________________________
string DCCSPPXSecWithCache::ProbeAsString (int probe_pdg, int probe_helicity) const
{
  string s_probe_pdg = std::to_string(probe_pdg);
  char c_hel = 0;
  if (probe_helicity == -1)
   c_hel = 'L';
  else if (probe_helicity == 1)
   c_hel = 'R';
  switch (probe_pdg)
  {
    case kPdgNuE:
    case kPdgNuMu:
    case kPdgNuTau:
        return ";v:" + s_probe_pdg;
    case kPdgAntiNuE:
    case kPdgAntiNuMu:
    case kPdgAntiNuTau:
        return ";vb:" + s_probe_pdg;
    case kPdgElectron:
    case kPdgPositron:
    case kPdgMuon:
    case kPdgAntiMuon:
    case kPdgTau:
    case kPdgAntiTau:
          return ";l:" + s_probe_pdg + c_hel;
  }
  return ";probe:" + s_probe_pdg + c_hel;
}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E::d3XSecSPP_dWQ2CosTheta_E(
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


  if (Enu < kps->Threshold())
  {
    isZero = true;
    return;
  }
  
  Wl  = kps->WLim_SPP();
  // model restrictions
  Wl.min  = TMath::Max (Wl.min,  1.08);
  Wl.max  = TMath::Max (Wl.max,  1.08);
  Wl.max  = TMath::Min (Wl.max,  2.00);

  if (fWcut >= Wl.min)
    Wl.max = TMath::Min(fWcut,Wl.max);
  
  
}
genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E::~d3XSecSPP_dWQ2CosTheta_E()
{

}
unsigned int genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E::DoEval(const double * xin) const
{

  // outputs:
  //   differential cross section [1/GeV^3] for Resonance single pion production production
  //

  if (isZero) return 0.;
  
  double W  = Wl.min + (Wl.max - Wl.min)*xin[0];
  fInteraction->KinePtr()->SetW(W);
   
  Range1D_t Q2l = kps->Q2Lim_W_SPP();
  // model restrictions
  Q2l.min = TMath::Max (Q2l.min, 0.00);
  Q2l.max = TMath::Min (Q2l.max, 3.00);
  
  
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[2]); //CosTheta
  
  
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpfE)*(Wl.max-Wl.min)*(Q2l.max-Q2l.min)*2;
  return xsec;
}
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSecSPP_dWQ2CosTheta_E(fModel,fInteraction,fWcut);
}
