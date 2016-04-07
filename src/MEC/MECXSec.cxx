//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "MEC/MECXSec.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"
#include "Utils/XSecSplineList.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
MECXSec::MECXSec() :
XSecIntegratorI("genie::MECXSec")
{

}
//____________________________________________________________________________
MECXSec::MECXSec(string config) :
XSecIntegratorI("genie::MECXSec", config)
{

}
//____________________________________________________________________________
MECXSec::~MECXSec()
{

}
//____________________________________________________________________________
double MECXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;
  
  const KPhaseSpace & kps = in->PhaseSpace(); // only OK phase space for this
  if(!kps.IsAboveThreshold()) {
     LOG("MECXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);
  
  // T, costh limits
  // NOTE: Where are these limits actually computed
  double Enu = in->InitState().ProbeE(kRfLab);
  double kine_min[2] = {  0,  -1 }; 
  double kine_max[2] = { Enu,  1 }; 

  double xsec = 0;

  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IBaseFunctionMultiDim * func = 
        new utils::gsl::d2Xsec_dTCosth(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);        
  ROOT::Math::IntegratorMultiDim ig(
    *func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
  xsec = ig.Integral(kine_min, kine_max); 

  delete func;
  delete interaction;   

  return xsec;
}
//____________________________________________________________________________
void MECXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  fGSLIntgType   = fConfig->GetStringDef("gsl-integration-type" ,    "vegas");
  fGSLMaxEval    = (unsigned int) fConfig->GetIntDef("gsl-max-evals", 20000);
  fGSLRelTol     = fConfig->GetDoubleDef("gsl-relative-tolerance",    0.01);
  fSplitIntegral = fConfig->GetBoolDef("split-integral",              true);
}
//_____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dTCosth::d2Xsec_dTCosth(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
//____________________________________________________________________________
genie::utils::gsl::d2Xsec_dTCosth::~d2Xsec_dTCosth()
{
  
}   
//____________________________________________________________________________
unsigned int genie::utils::gsl::d2Xsec_dTCosth::NDim(void) const
{
  return 2;
}
//____________________________________________________________________________
double genie::utils::gsl::d2Xsec_dTCosth::DoEval(const double * xin) const
{
// inputs:
//    T [GeV]
//    cos(theta)
// outputs:
//   differential cross section (hbar=c=1 units)
//
 
  double T     = xin[0];
  double costh = xin[1];

  Kinematics * kinematics = fInteraction->KinePtr();
  kinematics->SetKV(kKVTl, T);
  kinematics->SetKV(kKVctl, costh);
  
  double xsec = fModel->XSec(fInteraction, kPSTlctl); 
  return xsec;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2Xsec_dTCosth::Clone() const
{
  return
    new genie::utils::gsl::d2Xsec_dTCosth(fModel,fInteraction);
}
//____________________________________________________________________________


