//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/MECXSec.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Framework/Utils/XSecSplineList.h"

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
  double Enu = in->InitState().ProbeE(kRfLab);
  double LepMass = in->FSPrimLepton()->Mass();
  double TMax = Enu - LepMass;
  double TMin = 0.0;
  double CosthMax = 1.0;
  double CosthMin = -1.0;
  if (Enu < fQ3Max) {
    TMin = 0 ;
    CosthMin = -1 ; 
  } else {
    TMin = TMath::Sqrt(TMath::Power(LepMass, 2) + TMath::Power((Enu - fQ3Max), 2)) - LepMass;
    CosthMin = TMath::Sqrt(1 - TMath::Power((fQ3Max / Enu ), 2));
  }

  double kine_min[2] = { TMin,  CosthMin }; 
  double kine_max[2] = { TMax,  CosthMax }; 

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
  GetParam( "NSV-Q3Max", fQ3Max ) ;

  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("vegas") ) ;

  int max ;
  GetParamDef( "gsl-max-eval", max, 20000 ) ;
  fGSLMaxEval    = (unsigned int) max ;

  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 0.01 ) ;
  GetParamDef( "split-integral", fSplitIntegral, true ) ;

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


