//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/XSection/COHXSecAR.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
COHXSecAR::COHXSecAR() :
  XSecIntegratorI("genie::COHXSecAR"), 
  fHasPion( false ), fHasPhoton( false ), 
  fOmegaIntegral( false ), ftIntegral( false ), 
  fGammaLimits( nullptr ) 
{

}
//____________________________________________________________________________
COHXSecAR::COHXSecAR(string config) :
  XSecIntegratorI("genie::COHXSecAR", config), 
  fHasPion( false ), fHasPhoton( false ), 
  fOmegaIntegral( false ), ftIntegral( false ), 
  fGammaLimits( nullptr ) 
{

}
//____________________________________________________________________________
COHXSecAR::~COHXSecAR()
{

}

//____________________________________________________________________________

double COHXSecAR::Integrate( const XSecAlgorithmI * model, const Interaction * in) const {

  if      ( fHasPion )   return IntegratePion( model, in ) ;
  else if ( fHasPhoton ) return IntegratePhoton( model, in ) ;
  
  return 0. ;

}
//____________________________________________________________________________
double COHXSecAR::IntegratePion(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  const InitialState & init_state = in -> InitState();

  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSecAR", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  Range1D_t y_lim = kps.Limits(kKVy);

  // Check this
  double Enu      = init_state.ProbeE(kRfLab);
  double Elep_min = (1.-y_lim.max) * Enu;
  double Elep_max = (1.-y_lim.min) * Enu;

  LOG("COHXSecAR", pINFO)
       << "Lepton energy integration range = [" << Elep_min << ", " << Elep_max << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  double xsec = 0;
  if (fSplitIntegral) {
    utils::gsl::dXSec_dElep_AR_pion func(model, interaction, fGSLIntgType, fGSLRelTol, fGSLMaxEval);
    
    //~ ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;
    ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;

    double abstol = 1; // Pretty sure this parameter is unused by ROOT.
    int size = 1000;  // Max number of subintervals, won't reach nearly this.
    int rule = 2; // See https://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC283
                  // Rule 2 is 21 points min

    ROOT::Math::Integrator ig( func,ig_type,abstol,fGSLRelTol,size,rule);
    
    xsec = ig.Integral(Elep_min, Elep_max) * (1E-38 * units::cm2);
  }
  else {
    double zero    = kASmallNum;
    double pi      = kPi-kASmallNum ;
    double twopi   = 2*kPi-kASmallNum ;

    //~ ROOT::Math::IBaseFunctionMultiDim * func =
          //~ new utils::gsl::wrap::d5Xsec_dEldOmegaldOmegapi(model, interaction);
    //~ double kine_min[5] = { Elep_min, zero , zero    , zero, zero };
    //~ double kine_max[5] = { Elep_max, pi   , twopi   , pi  , twopi};

    ROOT::Math::IBaseFunctionMultiDim * func =
          new utils::gsl::d4Xsec_dEldThetaldOmegapi(model, interaction);
    double kine_min[4] = { Elep_min, zero , zero    , zero    };
    double kine_max[4] = { Elep_max, pi   , pi      , twopi   };

    ROOT::Math::IntegrationMultiDim::Type ig_type =
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);

    double abstol = 1; //We mostly care about relative tolerance.
    ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  }

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
double COHXSecAR::IntegratePhoton( const XSecAlgorithmI * model, const Interaction * in) const {

  //const InitialState & init_state = in -> InitState();

  if(! model->ValidProcess(in) ) { return 0.; }

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
    LOG("COHXSecAR", pDEBUG)  << "*** Below energy threshold";
    return 0;
  }
  
  // Check this
  Range1D_t e_gamma     = fGammaLimits -> EGamma( *in ) ;
  Range1D_t theta_gamma = fGammaLimits -> ThetaGamma( *in ) ;
  Range1D_t phi_gamma   = fGammaLimits -> PhiGamma( *in ) ;

  // LOG("COHXSecAR", pINFO)
  //      << "Lepton energy integration range = [" << Elep_min << ", " << Elep_max << "]";

  Interaction interaction(*in);
  interaction.SetBit(kISkipProcessChk);
  //interaction.SetBit(kISkipKinematicChk);
  
  
  //for the time begin the option of splitting the integral is not there for photon

  // if (fSplitIntegral) {
  //   utils::gsl::dXSec_dElep_AR * func =
  //     new utils::gsl::dXSec_dElep_AR(model, interaction, fGSLIntgType, fGSLRelTol, fGSLMaxEval);
    
  //   //~ ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;
  //   ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
    
  //   double abstol = 1; // Pretty sure this parameter is unused by ROOT.
  //   int size = 1000;  // Max number of subintervals, won't reach nearly this.
  //   int rule = 2; // See https://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC283
  //                 // Rule 2 is 21 points min
  //   ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,size,rule);
    
  //   xsec = ig.Integral(Elep_min, Elep_max) * (1E-38 * units::cm2);
  //   delete func;
  // }
  // else {

  ROOT::Math::IBaseFunctionMultiDim * func = nullptr ;
  double min_second, max_second ;
  if ( fOmegaIntegral ) { 
    func = new utils::gsl::d5Xsec_dEgdOmegaldOmegag(model, & interaction);
    Range1D_t theta_lep   = fGammaLimits -> ThetaLepton( *in ) ;
    min_second = cos( theta_lep.max ) ;
    max_second = cos( theta_lep.min ) ;
  }
  else if ( ftIntegral ) { 
    func =  new utils::gsl::d4Xsec_dEgdtdThetagdPhig(model, & interaction);
    Range1D_t t = fGammaLimits -> t( *in ) ;
    min_second = t.min ;
    max_second = t.max ;
  }
  else {
    func =  new utils::gsl::d4Xsec_dEgdThetaldThetagdPhig(model, & interaction); 
    Range1D_t theta_lep   = fGammaLimits -> ThetaLepton( *in ) ;
    min_second = theta_lep.min ;
    max_second = theta_lep.max ;
  }
    

  std::array<double,4> kine_min = { e_gamma.min, 
				    min_second, 
				    fOmegaIntegral ? cos( theta_gamma.max ) : theta_gamma.min, 
				    phi_gamma.min } ;
  std::array<double,4> kine_max = { e_gamma.max, 
				    max_second,
				    fOmegaIntegral ? cos( theta_gamma.min ) : theta_gamma.max, 
				    phi_gamma.max } ;
  
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  
  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
  double xsec = ig.Integral(kine_min.data(), kine_max.data()) * (1E-38 * units::cm2) ;
  
  if ( fOmegaIntegral ) xsec *= 2 * constants::kPi ;
  
  delete func;
  //  }
  
  return xsec;
}
//____________________________________________________________________________
void COHXSecAR::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("vegas") ) ;

  int max_eval;
  GetParamDef( "gsl-max-eval", max_eval, 4000 ) ;
  fGSLMaxEval    = (unsigned int) max_eval ;

  GetParamDef( "gsl-relative-tolerance", fGSLRelTol,  0.01) ;
  GetParamDef( "split-integral", fSplitIntegral, true ) ;

  GetParamDef( "IsCOHPion",  fHasPion,   false ) ;
  GetParamDef( "IsCOHGamma", fHasPhoton, false ) ;

  bool error = false ;

  if ( fHasPhoton ) {
    GetParam( "OmegaPhaseSpace", fOmegaIntegral ) ;

    GetParamDef( "tPhaseSpace", ftIntegral, false ) ;


    const Algorithm * temp = SubAlg( "IntegrationLimits" ) ;
    fGammaLimits = dynamic_cast<const COHGammaIntegrationLimits *>( temp ) ;
    if (! fGammaLimits ) {
      LOG( "COHXSecAR", pERROR ) << "Gamma integration limits subalgo failed to load" ;
      error = true ;
    }
  }

  if ( !fHasPion && !fHasPhoton ) {
    LOG( "COHXSecAR", pERROR ) << "No pion nor gamma option has been requested" ;
    error = true ;
  }

  if ( fHasPion && fHasPhoton ) {
    LOG( "COHXSecAR", pERROR ) << "Pion and Gamma options have been requested at the same time" ;
    error = true ;
  }

  if ( error ) {
    LOG( "COHXSecAR", pFATAL ) << "Invalid configuration. Exiting" ;
    exit( 78 ) ;
  } 

}
//_____________________________________________________________________________
