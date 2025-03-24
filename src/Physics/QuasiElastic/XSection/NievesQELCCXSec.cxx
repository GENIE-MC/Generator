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
#include "Physics/QuasiElastic/XSection/NievesQELCCPXSec.h"


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
  
  GetParamDef( "gsl-1dim-integration-type", fGSLIntgType, std::string("adaptive") );
  GetParamDef( "gsl-1dim-relative-tolerance", fGSLRelTol, 1e-5 );
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
  fInteraction(interaction),
  fRmax(Rmax)
{  
  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  fEnu = init_state.ProbeE(kRfHitNucRest);
  fml = interaction->FSPrimLepton()->Mass();
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
  // outputs:
  //   differential cross section [1/GeV^3] for Resonance single pion production production
  //
  double El   = fml + (fEnu - fml)*xin[0];
  double cost = -1 + 2*xin[1];
  double R    = fRmax*xin[2];
  
  double Pl   = TMath::Sqrt(El*El - fml*fml);
  double sint = TMath::Sqrt(1 - cost*cost);
  
  Kinematics * kinematics = fInteraction->KinePtr();
  kinematics->SetFSLeptonP4(Pl*sint, 0, Pl*cost, El);
  
  const NievesQELCCPXSec * xsec_model = dynamic_cast<const NievesQELCCPXSec*>(fModel);
  double xsec = xsec_model->IntegratedOverMomentum(fInteraction, R, 0);
  
  double J = (fEnu - fml)*2*fRmax;
  
  return xsec*J/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSec_dElepdCosThetalepdR_E(fModel,fInteraction,fRmax);
}
