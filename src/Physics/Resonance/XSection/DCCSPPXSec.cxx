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

#include <string>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

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
#include "Physics/Resonance/XSection/DCCSPPXSec.h"
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
DCCSPPXSec::DCCSPPXSec() :
DCCSPPXSecWithCache("genie::DCCSPPXSec")
{

}
//____________________________________________________________________________
DCCSPPXSec::DCCSPPXSec(string config) :
DCCSPPXSecWithCache("genie::DCCSPPXSec", config)
{

}
//____________________________________________________________________________
DCCSPPXSec::~DCCSPPXSec()
{

}
//____________________________________________________________________________
double DCCSPPXSec::Integrate(
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

  if (Enu < kps.Threshold()) return 0.;
  
  fSinglePionProductionXSecModel = model;

  InteractionType_t it    = proc_info.InteractionTypeId();
  int nucleon_pdgc        = target.HitNucPdg();
  int probe_pdgc          = init_state.ProbePdg();
  int probe_helicity      = init_state.ProbeHelicity();
  std::string nc_nuc   = this->ProbeAsString(probe_pdgc, probe_helicity);
  
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
      SLOG("DCCSPPXSec", pNOTICE)  
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
     string key = this->CacheBranchName(spp_channel, it, probe_pdgc, probe_helicity);
     LOG("DCCSPPXSec", pINFO) 
         << "Finding cache branch with key: " << key;
     CacheBranchFx * cache_branch =
         dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
     if(!cache_branch) {
        char c_hel = 0;
        if (probe_helicity == -1)
          c_hel = 'L';
        else if (probe_helicity == 1)
          c_hel = 'R';
        LOG("DCCSPPXSec", pWARN)  
           << "No cached ResSPP data for input probe with pdg = "
           << probe_pdgc << c_hel;
        LOG("DCCSPPXSec", pWARN)  
           << "Wait while computing/caching ResSPP production xsec first...";

        this->CacheResExcitationXSec(interaction); 

        LOG("DCCSPPXSec", pINFO) << "Done caching resonance xsec data";
        LOG("DCCSPPXSec", pINFO) 
               << "Finding newly created cache branch with key: " << key;
        cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
        assert(cache_branch);
     }
     const CacheBranchFx & cbranch = (*cache_branch);
  
    //-- Get cached resonance xsec
    //   (If E>Emax, assume xsec = xsec(Emax) - but do not evaluate the
    //    cross section spline at the end of its energy range)
    double rxsec = (Enu<fEMax-1) ? cbranch(Enu) : cbranch(fEMax-1);
               
    SLOG("DCCSPPXSec", pNOTICE)
      << "XSec[Channel: " << SppChannel::AsString(spp_channel) << nc_nuc
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
    LOG("DCCSPPXSec", pINFO)
          << "*** Integrating d^3 XSec/dWdQ^2dCosTheta for Ch: "
          << SppChannel::AsString(spp_channel) << " at Ev = " << Enu;
    
    ROOT::Math::IBaseFunctionMultiDim * func= new utils::gsl::d3XSecSPP_dWQ2CosTheta_E(model, interaction, fWcut);
    ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type,0,fGSLRelTol,fGSLMaxEval);
    ig.SetFunction(*func);
    double kine_min[3] = { 0., 0., 0.};
    double kine_max[3] = { 1., 1., 1.};
    double xsec = ig.Integral(kine_min, kine_max);
    delete func;
      
    SLOG("DCCSPPXSec", pNOTICE)
      << "XSec[Channel: " << SppChannel::AsString(spp_channel) << nc_nuc
      << "]  (E="<< Enu << " GeV) = " << xsec/(1E-38 *genie::units::cm2) << " x 1E-38 cm^2";
    
    return xsec;
  }
  return 0;
}
//____________________________________________________________________________
void DCCSPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPXSec::LoadConfig(void)
{
  
  bool good_conf = true ;

   // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 0.01 ) ;
  GetParamDef( "gsl-max-eval", fGSLMaxEval, 100000 ) ;
  GetParam("UsePauliBlockingForRES", fUsePauliBlocking);
  GetParamDef("Wcut", fWcut, -1.);
  // Get upper E limit on res xsec spline (=f(E)) before assuming xsec=const
  GetParamDef( "ESplineMax", fEMax, 500. ) ;

  if ( fEMax < 20. ) {

    LOG("DCCSPPXSec", pERROR) << "E max is required to be at least 20 GeV, you set " << fEMax << " GeV" ;
    good_conf = false ;
  }

  if ( ! good_conf ) {
    LOG("DCCSPPXSec", pFATAL)
      << "Invalid configuration: Exiting" ;
    
    // From the FreeBSD Library Functions Manual
    //
    // EX_CONFIG (78)   Something was found in an unconfigured or miscon-
    //                  figured state.
    
    exit( 78 ) ;
    
  }
  
}
//____________________________________________________________________________

