//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
 University of Sussex

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <limits>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/DarkNeutrino/XSection/COHDNuXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
COHDNuXSec::COHDNuXSec() :
XSecIntegratorI("genie::COHDNuXSec")
{

}
//____________________________________________________________________________
COHDNuXSec::COHDNuXSec(string config) :
XSecIntegratorI("genie::COHDNuXSec", config)
{

}
//____________________________________________________________________________
COHDNuXSec::~COHDNuXSec()
{

}
//____________________________________________________________________________
double COHDNuXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  if(!model->ValidProcess(in) ) return 0.;
  if(!in->PhaseSpace().IsAboveThreshold()) return 0.;

  Interaction interaction(*in);
  interaction.SetBit(kISkipProcessChk);
  interaction.SetBit(kISkipKinematicChk);

  ROOT::Math::IntegrationOneDim::Type ig_type =
    utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  
  utils::gsl::dXSec_dEDNu_E func( model, & interaction, fDNuMass );
  Range1D_t DNuEnergy = func.IntegrationRange();
  
  ROOT::Math::Integrator ig( func, ig_type, fGSLAbsTol, fGSLRelTol, fGSLMaxSizeOfSubintervals );
  double xsec = ig.Integral(DNuEnergy.min, DNuEnergy.max) * (1E-38 * units::cm2); // units: GeV^-2

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  LOG("COHDNuXSec", pINFO)
    << "XSec[COHDNu] (E = " << Ev << " GeV) = " << xsec/(units::cm2) << " cm^2";
  
  return xsec;
}
//____________________________________________________________________________
void COHDNuXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDNuXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDNuXSec::LoadConfig(void)
{
  fDNuMass = 0.;
  this->GetParam("Dark-NeutrinoMass", fDNuMass);

  // Get GSL integration type & relative tolerance
  this->GetParamDef("gsl-integration-type",   fGSLIntgType, string("adaptive"));
  this->GetParamDef("gsl-relative-tolerance", fGSLRelTol,   1E-2);
  this->GetParamDef("gsl-absolute-tolerance", fGSLAbsTol,   std::numeric_limits<double>::max());
  
  int max_subint = 0 ;
  this->GetParamDef("gsl-max-subintervals",   max_subint,   500);
  fGSLMaxSizeOfSubintervals = (unsigned int) max_subint ;

  // unused parameters from XSecIntegratorI
  fGSLMaxEval = fGSLMinEval = 0 ;
}
//____________________________________________________________________________
