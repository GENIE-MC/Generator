//_____________________________________________________________________________________
/*!

\namespace  genie::utils::gsl::wrap

\brief      GENIE differential cross section function wrappers for GSL integrators.

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            STFC, Rutherford Appleton Laboratory

\created    Sep 01, 2009

\cpright    Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//_____________________________________________________________________________________

#ifndef _GENIE_XSEC_FUNCTION_GSL_WRAPPERS_H_
#define _GENIE_XSEC_FUNCTION_GSL_WRAPPERS_H_

#include <Math/IFunction.h>

namespace genie {

class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {
namespace wrap  {

//.....................................................................................
//
// genie::utils::gsl::wrap::dXSec_dQ2_E
// A 1-D cross section function: dxsec/dQ2 = f(Q2)|(fixed E)
//
class dXSec_dQ2_E: public ROOT::Math::IBaseFunctionOneDim
{
public:
  dXSec_dQ2_E(const XSecAlgorithmI * m, const Interaction * i);
 ~dXSec_dQ2_E();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::dXSec_dy_E
// A 1-D cross section function: dxsec/dy = f(y)|(fixed E)
//
class dXSec_dy_E: public ROOT::Math::IBaseFunctionOneDim
{
public:
  dXSec_dy_E(const XSecAlgorithmI * m, const Interaction * i);
 ~dXSec_dy_E();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dxdy_E
// A 2-D cross section function: d2xsec/dxdy = f(x,y)|(fixed E)
//
class d2XSec_dxdy_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSec_dxdy_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d2XSec_dxdy_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dWdQ2_E
// A 2-D cross section function: d2xsec/dWdQ2 = f(W,Q2)|(fixed E)
//
class d2XSec_dWdQ2_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSec_dWdQ2_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d2XSec_dWdQ2_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dxdy_Ex
// A 1-D cross section function: d2xsec/dxdy = f(y)|(fixed:E,x)
//
class d2XSec_dxdy_Ex: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dxdy_Ex(const XSecAlgorithmI * m, const Interaction * i, double x);
 ~d2XSec_dxdy_Ex();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double fx;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dxdy_Ey 
// A 1-D cross section function: d2xsec/dxdy = f(x)|(fixed:E,y) 
//
class d2XSec_dxdy_Ey: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dxdy_Ey(const XSecAlgorithmI * m, const Interaction * i, double x);
 ~d2XSec_dxdy_Ey();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fy;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dWdQ2_EW
// A 1-D cross section function: d2xsec/dWdQ2= f(Q2)|(fixed:E,W)
//
class d2XSec_dWdQ2_EW: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dWdQ2_EW( const XSecAlgorithmI * m, const Interaction * i, double W);
 ~d2XSec_dWdQ2_EW();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fW;
};

//.....................................................................................
//
// genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2
// A 1-D cross section function: d2xsec/dWdQ2= f(W)|(fixed:E,Q2)
//
class d2XSec_dWdQ2_EQ2: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dWdQ2_EQ2(const XSecAlgorithmI * m, const Interaction * i, double Q2);
 ~d2XSec_dWdQ2_EQ2();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (const double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fQ2;
};

//.....................................................................................

} // wrap  namespace
} // gsl   namespace
} // utils namespace
} // genie namespace

#endif   

