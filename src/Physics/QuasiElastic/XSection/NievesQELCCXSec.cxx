//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin JINR
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/NievesQELCCXSec.h"

#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
NievesQELCCXSec::NievesQELCCXSec() : XSecIntegratorI("genie::NievesQELCCXSec")
{

}
//____________________________________________________________________________
NievesQELCCXSec::NievesQELCCXSec(std::string config) : XSecIntegratorI("genie::NievesQELCCXSec", config)
{

}
//____________________________________________________________________________
double NievesQELCCXSec::Integrate(const XSecAlgorithmI* model, const Interaction* in) const
{
    Target * tgt = in->InitStatePtr()->TgtPtr();
    double Rmax =  tgt->HitNucPosition();
    utils::gsl::d3XSec_dElepdCosThetalepdR_E func(model, in, Rmax ) ; 
    ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
    ROOT::Math::IntegratorMultiDim ig(ig_type, 0, fGSLRelTol, fGSLMaxEval);
    if (ig_type == ROOT::Math::IntegrationMultiDim::kADAPTIVE) 
    {
      ROOT::Math::AdaptiveIntegratorMultiDim * cast = dynamic_cast<ROOT::Math::AdaptiveIntegratorMultiDim*>( ig.GetIntegrator() );
      assert(cast);
      cast->SetMinPts(fGSLMinEval);
    }
    ig.SetFunction(func);
    double kine_min[3] = { 0., 0., 0.};
    double kine_max[3] = { 1., 1., 1.};
    double xsec = ig.Integral(kine_min, kine_max)*(1E-38 * units::cm2);
    return xsec;
}
//____________________________________________________________________________
void NievesQELCCXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NievesQELCCXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, std::string("adaptive") );
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-5 );
  int max;
  GetParamDef( "gsl-max-eval", max, 100000 );
  fGSLMaxEval  = static_cast<unsigned int>( max );
  int min;
  GetParamDef( "gsl-min-eval", min, 7500 ) ;
  fGSLMinEval  = static_cast<unsigned int>( min );
  
  GetParamDef( "gsl-1dim-integration-type", fGSL1DimIntgType, std::string("adaptive") );
  GetParamDef( "gsl-1dim-relative-tolerance", fGSL1DimRelTol, 1e-5 );
  GetParamDef( "gsl-1dim-max-eval", max, 100000 );
  fGSL1DimMaxEval  = static_cast<unsigned int>( max );
}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::d3XSec_dElepdCosThetalepdR_E(
       const XSecAlgorithmI * m, const Interaction * interaction, double  Rmax) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fXsec_model( dynamic_cast<const NievesQELCCPXSec*>(fModel) ),
  fInteraction(interaction),
  fRmax(Rmax)
{  
  // Get kinematical parameters
  AlgFactory * algf = AlgFactory::Instance();
  sm_utils = const_cast<genie::SmithMonizUtils *>(dynamic_cast<const genie::SmithMonizUtils *>(algf->GetAlgorithm("genie::SmithMonizUtils","Default")));
  const InitialState & init_state = interaction -> InitState();
  fEnu = init_state.ProbeE(kRfHitNucRest);
  fml = interaction->FSPrimLepton()->Mass();
  fml2 = fml*fml;
  sm_utils->SetInteraction(interaction);
  fKinematics = fInteraction->KinePtr();
}
genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::~d3XSec_dElepdCosThetalepdR_E()
{

}
unsigned int genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::DoEval(const double * xin) const
{
  // inputs:
  //    normalized Q2 from 0 to 1
  //    normalized v  from 0 to 1
  //    normalized R  from 0 to 1
  // outputs:
  //   differential cross section [10^-38 cm^2]
  //
  double R = fRmax*xin[2];
  double kFi, kFf;
  fXsec_model->ModelNuclParams(fInteraction, R, kFi, kFf);
  sm_utils->SetBindingEnergy(0);
  sm_utils->SetInitialFermiMomentum(kFi);
  sm_utils->SetFinalFermiMomentum(kFf);
  
  Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
  double Q2     = (rQ2.max - rQ2.min)*xin[0] + rQ2.min;
  Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
  double v      = (rv.max - rv.min)*xin[1] + rv.min;
  
  double El    = fEnu - v;
  if (El < fml) return 0.0;
  double Pl    = TMath::Sqrt(El*El - fml2);
  double cosTl = (El - (Q2 + fml2)/2/fEnu)/Pl;
  if (cosTl < -1.0 || cosTl > 1.0 ) return 0.0;
  double sinTl = TMath::Sqrt(1 - cosTl*cosTl);
  fKinematics->SetFSLeptonP4(Pl*sinTl, 0, Pl*cosTl, El);
  
  double xsec = fXsec_model->IntegratedOverMomentum(fInteraction, R, 0);
  
  // Jacobian for transformation d/dEldcosT->d/dQ2dv
  double J      = (rQ2.max - rQ2.min)*(rv.max - rv.min)*fRmax/2/fEnu/Pl; 
  
  return xsec*J/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E(fModel,fInteraction,fRmax);
}
