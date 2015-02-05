//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Chris Marshall and Martti Nirkko

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "AtharSingleKaon/ASKGSL.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
//#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/GSLUtils.h"

using namespace genie;

//____________________________________________________________________________
genie::utils::gsl::d3Xsec_dTldTkdCosThetal::d3Xsec_dTldTkdCosThetal(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
genie::utils::gsl::d3Xsec_dTldTkdCosThetal::~d3Xsec_dTldTkdCosThetal()
{
  
}   
unsigned int genie::utils::gsl::d3Xsec_dTldTkdCosThetal::NDim(void) const
{
  // phi_kq is important for kinematics generation
  // But dependence is weak so we will not use it in the integration
  return 3;
}
double genie::utils::gsl::d3Xsec_dTldTkdCosThetal::DoEval(const double * xin) const
{
// inputs:
//    Tl [GeV]
//    Tk [GeV]
//    cosine theta l
//    * calculate phi_kq based on neutrino energy -- this is for the integral only
// outputs:
//   differential cross section [10^-38 cm^2]
//
 
  double Enu = fInteraction->InitState().ProbeE(kRfLab);
    
  double phikq = 0.5*constants::kPi;
  if( Enu > 3.0 ) phikq = 0.55*constants::kPi;
  else if( Enu > 1.0 ) phikq = constants::kPi*(0.5 + 0.025*(Enu-1.0));

  Kinematics * kinematics = fInteraction->KinePtr();

  double T_l         = xin[0];
  double T_k         = xin[1];
  //double cos_theta_l = xin[2];
  double log_oneminuscostheta = xin[2];
  double cos_theta_l = 1.0 - TMath::Exp(log_oneminuscostheta);
  double J = 1.0 - cos_theta_l; // Jacobian for transformation
    
  kinematics->SetKV(kKVTl, T_l);
  kinematics->SetKV(kKVTk, T_k);
  kinematics->SetKV(kKVctl, cos_theta_l);
  kinematics->SetKV(kKVphikq, phikq);
  
  double xsec = fModel->XSec(fInteraction);
  LOG( "GXSecFunc", pDEBUG ) 
     << "t_l = " << T_l << " t_k = " << T_k << " costhetal = " << cos_theta_l << " phikq = " << phikq 
     << " enu = " << Enu << " Xsec = " << xsec;

  // integrate out phi_kq by multiplying by 2pi
  
  xsec *= 2.0 * genie::constants::kPi * J;
 
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d3Xsec_dTldTkdCosThetal::Clone() const
{
  return
    new genie::utils::gsl::d3Xsec_dTldTkdCosThetal(fModel,fInteraction);
}
//____________________________________________________________________________

