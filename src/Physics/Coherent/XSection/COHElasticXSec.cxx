//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author:  Costas Andreopoulos <costas.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/Coherent/XSection/COHElasticXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
COHElasticXSec::COHElasticXSec() :
XSecIntegratorI("genie::COHElasticXSec")
{

}
//____________________________________________________________________________
COHElasticXSec::COHElasticXSec(string config) :
XSecIntegratorI("genie::COHElasticXSec", config)
{

}
//____________________________________________________________________________
COHElasticXSec::~COHElasticXSec()
{

}
//____________________________________________________________________________
double COHElasticXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("CEvNS", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  double E = in->InitState().ProbeE(kRfLab);
  Range1D_t Q2;
  Q2.min = 0; // 1E-12; //controls::kASmallNum;
  Q2.max = 4*E*E;

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  LOG("CEvNS", pINFO)
      << "Q2 integration range = [" << Q2.min << ", " << Q2.max << "] GeV^2";

  ROOT::Math::IntegrationOneDim::Type ig_type =
          utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  ROOT::Math::IBaseFunctionOneDim * func =
      new utils::gsl::dXSec_dQ2_E(model, interaction);
  double abstol = 1; // We mostly care about relative tolerance
  ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(Q2.min, Q2.max) * (1E-38 * units::cm2); // units: GeV^-2

  delete func;

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  LOG("CEvNS", pINFO)
    << "XSec[COH] (E = " << Ev << " GeV) = " << xsec/(units::cm2) << " cm^2";

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void COHElasticXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElasticXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElasticXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  this->GetParamDef("gsl-integration-type",   fGSLIntgType, string("adaptive"));
  this->GetParamDef("gsl-relative-tolerance", fGSLRelTol,   1E-2);

  // max and minimum number of integrand evaluations
  int max_eval, min_eval ;
  this->GetParamDef("gsl-max-eval", max_eval, 500000);
  this->GetParamDef("gsl-min-eval", min_eval, 5000  );
  fGSLMaxEval = (unsigned int) max_eval ;
  fGSLMinEval = (unsigned int) min_eval ;

  // LOG("CEvNS", pDEBUG) << *this;
}
//____________________________________________________________________________
