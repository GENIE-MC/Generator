//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 04, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/GBuild.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "CrossSections/GXSecFunc.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
//____________________________________________________________________________
GXSecFunc::GXSecFunc(
         const XSecAlgorithmI * m, const Interaction * i, unsigned int ndim) :
GSFunc(ndim),
fModel(m),
fInteraction(i)
{
  assert (fModel);
  assert (fInteraction);
}
//____________________________________________________________________________
GXSecFunc::~GXSecFunc()
{
  fModel       = 0;
  fInteraction = 0;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DxDy_E::Integrand_D2XSec_DxDy_E(
                            const XSecAlgorithmI * m, const Interaction * i) :
GXSecFunc(m,i,2)
{

}
//____________________________________________________________________________
Integrand_D2XSec_DxDy_E::~Integrand_D2XSec_DxDy_E()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DxDy_E::operator() (const vector<double> & p)
{
  //-- set the kinematical peters
  double x = p[0];
  double y = p[1];
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(y);

  kinematics::UpdateWQ2FromXY(fInteraction);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_DXSec_DQ2_E::Integrand_DXSec_DQ2_E(
                      const XSecAlgorithmI * m, const Interaction * i) :
GXSecFunc(m,i,1)
{

}
//____________________________________________________________________________
Integrand_DXSec_DQ2_E::~Integrand_DXSec_DQ2_E()
{

}
//____________________________________________________________________________
double Integrand_DXSec_DQ2_E::operator() (const vector<double> & p)
{
  //-- set the kinematical peters
  double Q2 = p[0];
  fInteraction->KinePtr()->SetQ2(Q2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSQ2fE);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GXSecFunc", pDEBUG) << "xsec(Q2 = " << Q2 << ") = " << xsec;
#endif
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_E::Integrand_D2XSec_DWDQ2_E(
                           const XSecAlgorithmI * m, const Interaction * i) :
GXSecFunc(m,i,2)
{

}
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_E::~Integrand_D2XSec_DWDQ2_E()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DWDQ2_E::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double W  = p[0];
  double Q2 = p[1];
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

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_DXSec_Dy_E::Integrand_DXSec_Dy_E(
                           const XSecAlgorithmI * m, const Interaction * i) :
GXSecFunc(m,i,1)
{

}
//____________________________________________________________________________
Integrand_DXSec_Dy_E::~Integrand_DXSec_Dy_E()
{

}
//____________________________________________________________________________
double Integrand_DXSec_Dy_E::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double y = p[0];
  fInteraction->KinePtr()->Sety(y);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSyfE);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GXSecFunc", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
#endif
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DxDy_Ex::Integrand_D2XSec_DxDy_Ex(
                 const XSecAlgorithmI * m, const Interaction * i, double x) :
GXSecFunc(m,i,1)
{
  fx = x;
}
//____________________________________________________________________________
Integrand_D2XSec_DxDy_Ex::~Integrand_D2XSec_DxDy_Ex()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DxDy_Ex::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double y = p[0];

  fInteraction->KinePtr()->Setx(fx);
  fInteraction->KinePtr()->Sety(y);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DxDy_Ey::Integrand_D2XSec_DxDy_Ey(
                 const XSecAlgorithmI * m, const Interaction * i, double y) :
GXSecFunc(m,i,1)
{
  fy = y;
}
//____________________________________________________________________________
Integrand_D2XSec_DxDy_Ey::~Integrand_D2XSec_DxDy_Ey()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DxDy_Ey::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double x = p[0];

  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(fy);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_EW::Integrand_D2XSec_DWDQ2_EW(
                 const XSecAlgorithmI * m, const Interaction * i, double W) :
GXSecFunc(m,i,1)
{
  fW = W;
}
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_EW::~Integrand_D2XSec_DWDQ2_EW()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DWDQ2_EW::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double Q2 = p[0];

  fInteraction->KinePtr()->SetW(fW);
  fInteraction->KinePtr()->SetQ2(Q2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_EQ2::Integrand_D2XSec_DWDQ2_EQ2(
                 const XSecAlgorithmI * m, const Interaction * i, double Q2) :
GXSecFunc(m,i,1)
{
  fQ2 = Q2;
}
//____________________________________________________________________________
Integrand_D2XSec_DWDQ2_EQ2::~Integrand_D2XSec_DWDQ2_EQ2()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DWDQ2_EQ2::operator() (const vector<double> & p)
{
  //-- set the kinematical parameters
  double W = p[0];

  fInteraction->KinePtr()->SetW(W);
  fInteraction->KinePtr()->SetQ2(fQ2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction,kPSWQ2fE);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________


