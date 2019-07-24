//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          adapted from  fortran code provided by 
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, Joint Institute for Nuclear Research
          Vladimir Lyubushkin, Joint Institute for Nuclear Research
          Vadim Naumov <vnaumov@theor.jinr.ru>, Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
          
 For the class documentation see the corresponding header file.

 
*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/SmithMonizQELCCXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"

using namespace genie;
using std::ostringstream;

//____________________________________________________________________________
SmithMonizQELCCXSec::SmithMonizQELCCXSec() :
XSecIntegratorI("genie::SmithMonizQELCCXSec")
{

}
//____________________________________________________________________________
SmithMonizQELCCXSec::SmithMonizQELCCXSec(string config) :
XSecIntegratorI("genie::SmithMonizQELCCXSec", config)
{

}
//____________________________________________________________________________
SmithMonizQELCCXSec::~SmithMonizQELCCXSec()
{

}
//____________________________________________________________________________
double SmithMonizQELCCXSec::Integrate(
                  const XSecAlgorithmI * model, const Interaction * in) const
{
  LOG("SMQELXSec",pDEBUG) << "Beginning integrate";
  if(! model->ValidProcess(in)) return 0.;
  
  const InitialState & init_state = in -> InitState();
  const Target & target = init_state.Tgt();
  if (target.A()<3)
  {
     const KPhaseSpace & kps = in->PhaseSpace();
     if(!kps.IsAboveThreshold()) {
        LOG("SMQELXSec", pDEBUG)  << "*** Below energy threshold";
        return 0;
     }
     Range1D_t rQ2 = kps.Limits(kKVQ2);
     if(rQ2.min<0 || rQ2.max<0) return 0;
     Interaction * interaction = new Interaction(*in);
     interaction->SetBit(kISkipProcessChk);
     interaction->SetBit(kISkipKinematicChk);
     ROOT::Math::IBaseFunctionOneDim * func = new utils::gsl::dXSec_dQ2_E(model, interaction);
     ROOT::Math::IntegrationOneDim::Type ig_type = utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
     double abstol = 0; //We mostly care about relative tolerance
     ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxSizeOfSubintervals, fGSLRule);
     double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);
     delete func;
     delete interaction;
     return xsec;
  }
  else
  {   
     Interaction * interaction = new Interaction(*in);
     sm_utils->SetInteraction(in);
     if (interaction->InitState().ProbeE(kRfLab)<sm_utils->E_nu_thr_SM()) return 0;
     interaction->SetBit(kISkipProcessChk);
     interaction->SetBit(kISkipKinematicChk);
     double xsec = 0;
     
     
     ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d2Xsec_dQ2dv(model, interaction);
     double kine_min[2] = { 0, 0}; 
     double kine_max[2] = { 1, 1}; 
     
     ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType2D);
     
     double abstol = 0; //We mostly care about relative tolerance.
     ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol2D, fGSLMaxEval);
     
     xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
     delete func;
     delete interaction;
     
     return xsec;

  }
  
}
//____________________________________________________________________________
void SmithMonizQELCCXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizQELCCXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r("SmithMonizQELCCXSec_specific", false ) ;
  r.Set("sm_utils_algo", RgAlg("genie::SmithMonizUtils","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizQELCCXSec::LoadConfig(void)
{
  
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("gauss") );
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-3 );
  int max_size_of_subintervals;
  GetParamDef( "gsl-max-size-of-subintervals", max_size_of_subintervals, 40000);
  fGSLMaxSizeOfSubintervals = (unsigned int) max_size_of_subintervals;
  int rule;
  GetParamDef( "gsl-rule", rule, 3);
  fGSLRule = (unsigned int) rule;
  if (fGSLRule>6) fGSLRule=3;
  GetParamDef( "gsl-integration-type-2D", fGSLIntgType2D, string("adaptive") );
  GetParamDef( "gsl-relative-tolerance-2D", fGSLRelTol2D, 1e-7);
  GetParamDef( "gsl-max-eval", fGSLMaxEval, 1000000000);
  
  sm_utils = const_cast<genie::SmithMonizUtils *>(
               dynamic_cast<const genie::SmithMonizUtils *>(
                 this -> SubAlg("sm_utils_algo") ) );
}


//_____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dQ2dv::d2Xsec_dQ2dv(const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(interaction)
{
  AlgFactory * algf = AlgFactory::Instance();
  sm_utils = const_cast<genie::SmithMonizUtils *>(dynamic_cast<const genie::SmithMonizUtils *>(algf->GetAlgorithm("genie::SmithMonizUtils","Default")));
  sm_utils->SetInteraction(interaction);
}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dQ2dv::~d2Xsec_dQ2dv()
{
  
}   
//____________________________________________________________________________
unsigned int genie::utils::gsl::d2Xsec_dQ2dv::NDim(void) const
{
  return 2;
}
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dQ2dv::DoEval(const double * xin) const
{
// inputs:
//    normalized Q2 from 0 to 1
//    normalized v  from 0 to 1
// outputs:
//   differential cross section [10^-38 cm^2]
//
 
  Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
  double Q2     = (rQ2.max-rQ2.min)*xin[0]+rQ2.min;
  Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
  double v      = (rv.max-rv.min)*xin[1]+rv.min;
  double J  = (rQ2.max-rQ2.min)*(rv.max-rv.min); // Jacobian for transformation
    
  Kinematics * kinematics = fInteraction->KinePtr();
  kinematics->SetKV(kKVQ2, Q2);
  kinematics->SetKV(kKVv, v);
   
  double xsec=fModel->XSec(fInteraction, kPSQ2vfE);
  
  xsec *= J;
   
  return xsec/(1E-38 * units::cm2);
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2Xsec_dQ2dv::Clone() const
{
  return new genie::utils::gsl::d2Xsec_dQ2dv(fModel, fInteraction);
}


