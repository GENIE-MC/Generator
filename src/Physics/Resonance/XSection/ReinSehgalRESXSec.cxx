//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
 
 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.
 @ Jan 29, 2013 - CA
   Don't look-up depreciated $GDISABLECACHING environmental variable.
   Use the RunOpt singleton instead.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Physics/Resonance/XSection/ReinSehgalRESXSec.h"

using namespace genie;
using namespace genie::constants;
//using namespace genie::units;

//____________________________________________________________________________
ReinSehgalRESXSec::ReinSehgalRESXSec() :
ReinSehgalRESXSecWithCache("genie::ReinSehgalRESXSec")
{

}
//____________________________________________________________________________
ReinSehgalRESXSec::ReinSehgalRESXSec(string config) :
ReinSehgalRESXSecWithCache("genie::ReinSehgalRESXSec", config)
{

}
//____________________________________________________________________________
ReinSehgalRESXSec::~ReinSehgalRESXSec()
{

}
//____________________________________________________________________________
double ReinSehgalRESXSec::Integrate(
          const XSecAlgorithmI * model, const Interaction * interaction) const
{
  if(! model->ValidProcess(interaction) ) return 0.;
  fSingleResXSecModel = model;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("ReinSehgalRESXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  //-- Get init state and process information
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  InteractionType_t it = proc_info.InteractionTypeId();
  int nucleon_pdgc = target.HitNucPdg();
  int nu_pdgc      = init_state.ProbePdg();

  //-- Get neutrino energy in the struck nucleon rest frame
  double Ev = init_state.ProbeE(kRfHitNucRest);

  //-- Get the requested resonance
  Resonance_t res = interaction->ExclTag().Resonance();

  // If the input interaction is off a nuclear target, then chek whether
  // the corresponding free nucleon cross section already exists at the
  // cross section spline list.
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) {
    Interaction * in = new Interaction(*interaction);
    if(pdg::IsProton(nucleon_pdgc)) { 
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeP); 
    } else { 
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeN); 
    }
    if(xsl->SplineExists(model,in)) {
      const Spline * spl = xsl->GetSpline(model, in);
      double xsec = spl->Evaluate(Ev);
      SLOG("ReinSehgalResT", pNOTICE)  
         << "XSec[RES/" << utils::res::AsString(res)<< "/free] (Ev = " 
         << Ev << " GeV) = " << xsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";
      if(! interaction->TestBit(kIAssumeFreeNucleon) ) {
        int NNucl = (pdg::IsProton(nucleon_pdgc)) ? target.Z() : target.N();
        xsec *= NNucl;
      }
      delete in;
      return xsec;
    }
    delete in;
  }

  // There was no corresponding free nucleon spline saved in XSecSplineList that
  // could be used to speed up this calculation.
  // Check whether local caching of free nucleon cross sections is allowed.
  // If yes, store free nucleon cross sections at a cache branch and use those
  // at any subsequent call.
  //
  bool bare_xsec_pre_calc = RunOpt::Instance()->BareXSecPreCalc();
  if(bare_xsec_pre_calc && !fUsePauliBlocking) {
     Cache * cache = Cache::Instance();
     string key = this->CacheBranchName(res, it, nu_pdgc, nucleon_pdgc);
     LOG("ReinSehgalResT", pINFO) 
         << "Finding cache branch with key: " << key;
     CacheBranchFx * cache_branch =
         dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
     if(!cache_branch) {
        LOG("ReinSehgalResT", pWARN)  
           << "No cached RES v-production data for input neutrino"
           << " (pdgc: " << nu_pdgc << ")";
        LOG("ReinSehgalResT", pWARN)  
           << "Wait while computing/caching RES production xsec first...";

        this->CacheResExcitationXSec(interaction); 

        LOG("ReinSehgalResT", pINFO) << "Done caching resonance xsec data";
        LOG("ReinSehgalResT", pINFO) 
               << "Finding newly created cache branch with key: " << key;
        cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
        assert(cache_branch);
     }
     const CacheBranchFx & cbranch = (*cache_branch);
  
    //-- Get cached resonance neutrinoproduction xsec
    //   (If E>Emax, assume xsec = xsec(Emax) - but do not evaluate the
    //    cross section spline at the end of its energy range-)
    double rxsec = (Ev<fEMax-1) ? cbranch(Ev) : cbranch(fEMax-1);

    SLOG("ReinSehgalResT", pNOTICE)  
       << "XSec[RES/" << utils::res::AsString(res)<< "/free] (Ev = " 
       << Ev << " GeV) = " << rxsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";

     if( interaction->TestBit(kIAssumeFreeNucleon) ) return rxsec;

     int NNucl = (pdg::IsProton(nucleon_pdgc)) ? target.Z() : target.N();
     rxsec*=NNucl; // nuclear xsec 
     return rxsec;
  } // disable local caching

  // Just go ahead and integrate the input differential cross section for the
  // specified interaction.  
  else {

    Range1D_t rW  = kps.Limits(kKVW);
    Range1D_t rQ2 = kps.Limits(kKVQ2);

    LOG("ReinSehgalResC", pINFO)
          << "*** Integrating d^2 XSec/dWdQ^2 for R: "
          << utils::res::AsString(res) << " at Ev = " << Ev;
    LOG("ReinSehgalResC", pINFO)
          << "{W}   = " << rW.min  << ", " << rW.max;
    LOG("ReinSehgalResC", pINFO)
          << "{Q^2} = " << rQ2.min << ", " << rQ2.max;

    ROOT::Math::IBaseFunctionMultiDim * func =
        new utils::gsl::d2XSec_dWdQ2_E(model, interaction);
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
        utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval); 
    ig.SetFunction(*func);
    double kine_min[2] = { rW.min, rQ2.min };
    double kine_max[2] = { rW.max, rQ2.max };
    double xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);

    delete func;
    return xsec;
  }
  return 0;
}
//____________________________________________________________________________
void ReinSehgalRESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalRESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalRESXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 0.01 ) ;
  GetParamDef( "gsl-max-eval", fGSLMaxEval, 100000 ) ;
  GetParam("UsePauliBlockingForRES", fUsePauliBlocking);
  // Get upper E limit on res xsec spline (=f(E)) before assuming xsec=const
  GetParamDef( "ESplineMax", fEMax, 100. ) ;
  fEMax = TMath::Max(fEMax, 20.); // don't accept user Emax if less than 20 GeV

  // Create the baryon resonance list specified in the config.
  fResList.Clear();
  string resonances ;
  GetParam( "ResonanceNameList", resonances ) ;
  fResList.DecodeFromNameList(resonances);

}
//____________________________________________________________________________
