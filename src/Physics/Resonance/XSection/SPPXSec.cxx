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

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/AdaptiveIntegratorMultiDim.h>

#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/Resonance/XSection/SPPXSec.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Numerical/GSLUtils.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::units;

//____________________________________________________________________________
SPPXSec::SPPXSec() :
SPPXSecWithCache("genie::SPPXSec")
{

}
//____________________________________________________________________________
SPPXSec::SPPXSec(string config) :
SPPXSecWithCache("genie::SPPXSec", config)
{

}
//____________________________________________________________________________
SPPXSec::~SPPXSec()
{

}
//____________________________________________________________________________
double SPPXSec::Integrate(
          const XSecAlgorithmI * model, const Interaction * interaction) const
{
  if(!model->ValidProcess(interaction) ) return 0.;

  //-- Get init state and process information
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();
  
  const KPhaseSpace& kps = interaction->PhaseSpace();
  
  double Enu = init_state.ProbeE(kRfHitNucRest);
  
  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  if (Enu < kps.Threshold_SPP_iso()) return 0.;
  
  fSinglePionProductionXSecModel = model;

  InteractionType_t it = proc_info.InteractionTypeId();
  int nucleon_pdgc = target.HitNucPdg();
  int nu_pdgc      = init_state.ProbePdg();
  string nc_nuc   = ((pdg::IsNeutrino(nu_pdgc)) ? ";v:" : ";vb:");
  
  // If the input interaction is off a nuclear target, then chek whether
  // the corresponding free nucleon cross section already exists at the
  // cross section spline list.
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) 
  {
    Interaction * in = new Interaction(*interaction);
    if(pdg::IsProton(nucleon_pdgc))  
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeP); 
    else 
      in->InitStatePtr()->TgtPtr()->SetId(kPdgTgtFreeN); 
    
    if(xsl->SplineExists(model,in)) 
    {
      const Spline * spl = xsl->GetSpline(model, in);
      double xsec = spl->Evaluate(Enu);
      SLOG("SPPXSec", pNOTICE)  
         << "XSec[Channel/" << SppChannel::AsString(spp_channel) << "/free] (Ev = " 
               << Enu << " GeV) = " << xsec/(1E-38 *cm2)<< " x 1E-38 cm^2";
      if(! interaction->TestBit(kIAssumeFreeNucleon) ) 
      {
        int NNucl = (SppChannel::InitStateNucleon(spp_channel) == kPdgProton) ? target.Z() : target.N();
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
  if(bare_xsec_pre_calc && !fUsePauliBlocking) 
  {
     Cache * cache = Cache::Instance();
     string key = this->CacheBranchName(spp_channel, it, nu_pdgc);
     LOG("SPPXSec", pINFO) 
         << "Finding cache branch with key: " << key;
     CacheBranchFx * cache_branch =
         dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
     if(!cache_branch) {
        LOG("SPPXSec", pWARN)  
           << "No cached ResSPP v-production data for input neutrino"
           << " (pdgc: " << nu_pdgc << ")";
        LOG("SPPXSec", pWARN)  
           << "Wait while computing/caching ResSPP production xsec first...";

        this->CacheResExcitationXSec(interaction); 

        LOG("SPPXSec", pINFO) << "Done caching resonance xsec data";
        LOG("SPPXSec", pINFO) 
               << "Finding newly created cache branch with key: " << key;
        cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
        assert(cache_branch);
     }
     const CacheBranchFx & cbranch = (*cache_branch);
  
    //-- Get cached resonance neutrinoproduction xsec
    //   (If E>Emax, assume xsec = xsec(Emax) - but do not evaluate the
    //    cross section spline at the end of its energy range-)
    double rxsec = (Enu<fEMax-1) ? cbranch(Enu) : cbranch(fEMax-1);

               
    SLOG("SPPXSec", pNOTICE)
      << "XSec[Channel: " << SppChannel::AsString(spp_channel) << nc_nuc  << nu_pdgc
      << "/free]  (E="<< Enu << " GeV) = " << rxsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";

     if( interaction->TestBit(kIAssumeFreeNucleon) ) return rxsec;

     int NNucl = (pdg::IsProton(nucleon_pdgc)) ? target.Z() : target.N();
     rxsec*=NNucl; // nuclear xsec 
     return rxsec;
  } // disable local caching

  // Just go ahead and integrate the input differential cross section for the
  // specified interaction.  
  else 
  {
    LOG("SPPXSec", pINFO)
          << "*** Integrating d^3 XSec/dWdQ^2dCosTheta for Ch: "
          << SppChannel::AsString(spp_channel) << " at Ev = " << Enu;
    
    ROOT::Math::IBaseFunctionMultiDim * func= new utils::gsl::d3XSecMK_dWQ2CosTheta_E(model, interaction, fWcut);
    ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
    if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) 
    {
      ROOT::Math::AdaptiveIntegratorMultiDim * cast = dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
      assert(cast);
      cast->SetMinPts(fGSLMinEval);
    }
    ig.SetFunction(*func);
    double kine_min[3] = { 0., 0., 0.};
    double kine_max[3] = { 1., 1., 1.};
    double xsec = ig.Integral(kine_min, kine_max)*(1E-38 * units::cm2);
    delete func;
      
    SLOG("SPPXSec", pNOTICE)
      << "XSec[Channel: " << SppChannel::AsString(spp_channel) << nc_nuc  << nu_pdgc
      << "]  (E="<< Enu << " GeV) = " << xsec << " x 1E-38 cm^2";
    return xsec;
  }
  return 0;
}
//____________________________________________________________________________
void SPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SPPXSec::LoadConfig(void)
{
  
  bool good_conf = true ;

   // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-5 ) ;
  GetParamDef( "gsl-max-eval", fGSLMaxEval, 1000000 );
  GetParamDef( "gsl-min-eval", fGSLMinEval, 7500 ) ;
  GetParam("UsePauliBlockingForRES", fUsePauliBlocking);
  GetParamDef("Wcut", fWcut, -1.);
  // Get upper Emax limit on cached free nucleon xsec spline, 
  // after this value it assume that xsec=xsec(Emax)
  GetParamDef( "ESplineMax", fEMax, 250. ) ;

  if ( fEMax < 20. ) {

    LOG("SPPXSec", pERROR) << "E max is required to be at least 20 GeV, you set " << fEMax << " GeV" ;
    good_conf = false ;
  }

  if ( ! good_conf ) {
    LOG("SPPXSec", pFATAL)
      << "Invalid configuration: Exiting" ;
    
    // From the FreeBSD Library Functions Manual
    //
    // EX_CONFIG (78)   Something was found in an unconfigured or miscon-
    //                  figured state.
    
    exit( 78 ) ;
    
  }
  
}
//____________________________________________________________________________

