//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 03, 2009 - CA
   Was first added in v2.5.1

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;

//____________________________________________________________________________
genie::utils::gsl::wrap::dXSec_dQ2_E::dXSec_dQ2_E(
    const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::wrap::dXSec_dQ2_E::~dXSec_dQ2_E()
{

}
unsigned int genie::utils::gsl::wrap::dXSec_dQ2_E::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::dXSec_dQ2_E::DoEval(const double xin) const
{
// inputs:  
//    Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2 / GeV^2]
//
  double Q2 = xin;
  fInteraction->KinePtr()->SetQ2(Q2);
  double xsec = fModel->XSec(fInteraction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GXSecFunc", pDEBUG) << "xsec(Q2 = " << Q2 << ") = " << xsec;
#endif
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::dXSec_dQ2_E::Clone() const
{
  return
    new genie::utils::gsl::wrap::dXSec_dQ2_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::dXSec_dy_E::dXSec_dy_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::wrap::dXSec_dy_E::~dXSec_dy_E()
{

}
unsigned int genie::utils::gsl::wrap::dXSec_dy_E::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::dXSec_dy_E::DoEval(const double xin) const
{
// inputs:  
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double y = xin;
  fInteraction->KinePtr()->Sety(y);
  double xsec = fModel->XSec(fInteraction, kPSyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GXSecFunc", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
#endif
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::dXSec_dy_E::Clone() const
{
  return
    new genie::utils::gsl::wrap::dXSec_dy_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dxdy_E::d2XSec_dxdy_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
genie::utils::gsl::wrap::d2XSec_dxdy_E::~d2XSec_dxdy_E()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dxdy_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::wrap::d2XSec_dxdy_E::DoEval(const double * xin) const
{
// inputs:  
//    x [-]
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double x = xin[0];
  double y = xin[1];
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::utils::gsl::wrap::d2XSec_dxdy_E::Clone() const
{
  return 
    new genie::utils::gsl::wrap::d2XSec_dxdy_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dWdQ2_E::d2XSec_dWdQ2_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::wrap::d2XSec_dWdQ2_E::~d2XSec_dWdQ2_E()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dWdQ2_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::wrap::d2XSec_dWdQ2_E::DoEval(const double * xin) const
{
// inputs:  
//    W  [GeV]
//    Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double W  = xin[0];
  double Q2 = xin[1];
  fInteraction->KinePtr()->SetW(W);
  fInteraction->KinePtr()->SetQ2(Q2);
  if(fInteraction->ProcInfo().IsDeepInelastic()) {
    double x=0,y=0;
    double E = fInteraction->InitState().ProbeE(kRfHitNucRest);
    double M = fInteraction->InitState().Tgt().HitNucP4Ptr()->M();
    kinematics::WQ2toXY(E,M,W,Q2,x,y);
    fInteraction->KinePtr()->Setx(x);
    fInteraction->KinePtr()->Sety(y);
  }
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::wrap::d2XSec_dWdQ2_E::Clone() const
{
  return 
    new genie::utils::gsl::wrap::d2XSec_dWdQ2_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dxdy_Ex::d2XSec_dxdy_Ex(
     const XSecAlgorithmI * m, const Interaction * i, double x) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fx(x)
{

}
genie::utils::gsl::wrap::d2XSec_dxdy_Ex::~d2XSec_dxdy_Ex()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dxdy_Ex::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::d2XSec_dxdy_Ex::DoEval(const double xin) const
{
// inputs:  
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double y = xin;
  fInteraction->KinePtr()->Setx(fx);
  fInteraction->KinePtr()->Sety(y);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::d2XSec_dxdy_Ex::Clone() const
{
  return
    new genie::utils::gsl::wrap::d2XSec_dxdy_Ex(fModel,fInteraction,fx);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dxdy_Ey::d2XSec_dxdy_Ey(
     const XSecAlgorithmI * m, const Interaction * i, double y) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fy(y)
{

}
genie::utils::gsl::wrap::d2XSec_dxdy_Ey::~d2XSec_dxdy_Ey()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dxdy_Ey::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::d2XSec_dxdy_Ey::DoEval(const double xin) const
{
// inputs:  
//    x [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double x = xin;
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(fy);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::d2XSec_dxdy_Ey::Clone() const
{
  return
    new genie::utils::gsl::wrap::d2XSec_dxdy_Ey(fModel,fInteraction,fy);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dWdQ2_EW::d2XSec_dWdQ2_EW(
     const XSecAlgorithmI * m, const Interaction * i, double W) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fW(W)
{

}
genie::utils::gsl::wrap::d2XSec_dWdQ2_EW::~d2XSec_dWdQ2_EW()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dWdQ2_EW::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::d2XSec_dWdQ2_EW::DoEval(const double xin) const
{
// inputs:  
//   Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double Q2 = xin;
  fInteraction->KinePtr()->SetW(fW);
  fInteraction->KinePtr()->SetQ2(Q2);
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::d2XSec_dWdQ2_EW::Clone(void) const
{
  return
    new genie::utils::gsl::wrap::d2XSec_dWdQ2_EW(fModel,fInteraction,fW);
}
//____________________________________________________________________________
genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2::d2XSec_dWdQ2_EQ2(
     const XSecAlgorithmI * m, const Interaction * i, double Q2) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fQ2(Q2)
{

}
genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2::~d2XSec_dWdQ2_EQ2()
{

}
unsigned int genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2::DoEval(const double xin) const
{
// inputs:  
//   W [GeV]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double W = xin;
  fInteraction->KinePtr()->SetW(W);
  fInteraction->KinePtr()->SetQ2(fQ2);
  double xsec = fModel->XSec(fInteraction,kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2::Clone() const
{
  return 
   new genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2(fModel,fInteraction,fQ2);
}
//____________________________________________________________________________
