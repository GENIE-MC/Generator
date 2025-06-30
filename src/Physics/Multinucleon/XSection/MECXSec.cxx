//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

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
#include "Physics/Multinucleon/XSection/MECUtils.h"
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

  Interaction interaction(*in);
  interaction.SetBit(kISkipProcessChk);
  interaction.SetBit(kISkipKinematicChk);

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
  genie::utils::mec::gsl::d2Xsec_dTCosth func(model, interaction, Enu, LepMass );
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  ROOT::Math::IntegratorMultiDim ig(func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  xsec = ig.Integral(kine_min, kine_max);

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
