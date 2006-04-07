#include <cassert>

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "CrossSections/GXSecFunc.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;

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
  fInteraction->GetKinematicsPtr()->Setx(x);
  fInteraction->GetKinematicsPtr()->Sety(y);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________
Integrand_D2XSec_DxDy_E_WQ2Cuts::Integrand_D2XSec_DxDy_E_WQ2Cuts(
            const XSecAlgorithmI * m, const Interaction * i,
                                         Range1D_t WCuts, Range1D_t Q2Cuts) :
GXSecFunc(m,i,2)
{
  fWCuts.Copy(WCuts);
  fQ2Cuts.Copy (Q2Cuts);
}
//____________________________________________________________________________
Integrand_D2XSec_DxDy_E_WQ2Cuts::~Integrand_D2XSec_DxDy_E_WQ2Cuts()
{

}
//____________________________________________________________________________
double Integrand_D2XSec_DxDy_E_WQ2Cuts::operator() (const vector<double> & p)
{
  //-- set the kinematical peters
  double x = p[0];
  double y = p[1];
  fInteraction->GetKinematicsPtr()->Setx(x);
  fInteraction->GetKinematicsPtr()->Sety(y);

  // get physical integration range for W and Q2
  Range1D_t rW  = utils::kinematics::WRange     (fInteraction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_xy (fInteraction);

  LOG("GXSecFunc", pDEBUG)
       << "\n Physical integration range: "
             << "W = [" << rW.min << ", " << rW.max << "] GeV "
                  << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // apply cuts
  utils::kinematics::ApplyCutsToKineLimits(rW,  fWCuts.min,   fWCuts.max );
  utils::kinematics::ApplyCutsToKineLimits(rQ2, fQ2Cuts.min,  fQ2Cuts.max);

  LOG("GXSecFunc", pDEBUG)
       << "\n Physical && User integration range: "
            << "W = [" << rW.min << ", " << rW.max << "] GeV "
                 << "Q2 = [" << rQ2.min << ", " << rQ2.max << "] GeV^2";

  // current W, Q2
  const InitialState & init_state = fInteraction -> GetInitialState();

  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonP4()->M();
  double M2 = TMath::Power(M,2);

  double currW2 = TMath::Max(0., M2 + 2*Ev*M*y*(1-x));
  double currW  = TMath::Sqrt(currW2);
  double currQ2 = TMath::Max(0., 2*M*Ev*x*y);

  bool in_range = utils::math::IsWithinLimits(currQ2, rQ2)
                             && utils::math::IsWithinLimits(currW, rW);
  if(!in_range) {
      LOG("GXSecFunc", pDEBUG) << "*** excluding from phase space: ";
      return 0;
  }

  LOG("GXSecFunc", pDEBUG) << "*** including in phase space: ";
  LOG("GXSecFunc", pDEBUG)
          << "x = " << x << ", y = " << y << ", Q2 = " << currQ2
                                << ", W = " << currW << ", Ev = " << Ev;
  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
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
  fInteraction->GetKinematicsPtr()->SetQ2(Q2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);

  LOG("GXSecFunc", pDEBUG) << "xsec(Q2 = " << Q2 << ") = " << xsec;
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
  //-- set the kinematical peters
  double W  = p[0];
  double Q2 = p[1];
  fInteraction->GetKinematicsPtr()->SetW(W);
  fInteraction->GetKinematicsPtr()->SetQ2(Q2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
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
  //-- set the kinematical peters
  double y = p[0];
  fInteraction->GetKinematicsPtr()->Sety(y);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);

  LOG("GXSecFunc", pDEBUG) << "xsec(y = " << y << ") = " << xsec;

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
  //-- set the kinematical peters
  double y = p[0];

  fInteraction->GetKinematicsPtr()->Setx(fx);
  fInteraction->GetKinematicsPtr()->Sety(y);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
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
  //-- set the kinematical peters
  double x = p[0];

  fInteraction->GetKinematicsPtr()->Setx(x);
  fInteraction->GetKinematicsPtr()->Sety(fy);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
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
  //-- set the kinematical peters
  double Q2 = p[0];

  fInteraction->GetKinematicsPtr()->SetW(fW);
  fInteraction->GetKinematicsPtr()->SetQ2(Q2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
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
  //-- set the kinematical peters
  double W = p[0];

  fInteraction->GetKinematicsPtr()->SetW(W);
  fInteraction->GetKinematicsPtr()->SetQ2(fQ2);

  //-- compute the cross section
  double xsec = fModel->XSec(fInteraction);
  return xsec;
}
//____________________________________________________________________________
//____________________________________________________________________________


