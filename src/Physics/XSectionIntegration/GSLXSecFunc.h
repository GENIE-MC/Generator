//_____________________________________________________________________________________
/*!

\namespace  genie::utils::gsl

\brief      GENIE differential cross section function wrappers for GSL integrators.

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    Sep 01, 2009

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//_____________________________________________________________________________________

#ifndef _GENIE_XSEC_FUNCTION_GSL_WRAPPERS_H_
#define _GENIE_XSEC_FUNCTION_GSL_WRAPPERS_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

namespace genie {

class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {

//.....................................................................................
//
// genie::utils::gsl::dXSec_dQ2_E
// A 1-D cross section function: dxsec/dQ2 = f(Q2)|(fixed E)
//
class dXSec_dQ2_E: public ROOT::Math::IBaseFunctionOneDim
{
public:
  dXSec_dQ2_E(const XSecAlgorithmI * m, const Interaction * i);
 ~dXSec_dQ2_E();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::dXSec_dy_E
// A 1-D cross section function: dxsec/dy = f(y)|(fixed E)
//
class dXSec_dy_E: public ROOT::Math::IBaseFunctionOneDim
{
public:
  dXSec_dy_E(const XSecAlgorithmI * m, const Interaction * i);
 ~dXSec_dy_E();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

//.....................................................................................
//
// genie::utils::gsl::d2XSec_dxdy_E
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
// genie::utils::gsl::d2XSec_dQ2dy_E
// A 2-D cross section function: d2xsec/dQ2dy = f(Q^2,y)|(fixed E)
//
class d2XSec_dQ2dy_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSec_dQ2dy_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d2XSec_dQ2dy_E();

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
// genie::utils::gsl::d2XSec_dQ2dydt_E
// A 3-D cross section function: d3xsec/dQ2dydt = f(Q^2,y,t)|(fixed E)
//
class d2XSec_dQ2dydt_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d2XSec_dQ2dydt_E(const XSecAlgorithmI * m, const Interaction * i);
 ~d2XSec_dQ2dydt_E();

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
// genie::utils::gsl::d3XSec_dxdydt_E
// A 3-D cross section function: d3xsec/dxdydt = f(x,y,t)|(fixed E)
//
class d3XSec_dxdydt_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3XSec_dxdydt_E(const XSecAlgorithmI * m, const Interaction * i);
  ~d3XSec_dxdydt_E();

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
// genie::utils::gsl::d2XSec_dWdQ2_E
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
// genie::utils::gsl::d2XSec_dxdy_Ex
// A 1-D cross section function: d2xsec/dxdy = f(y)|(fixed:E,x)
//
class d2XSec_dxdy_Ex: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dxdy_Ex(const XSecAlgorithmI * m, const Interaction * i, double x);
 ~d2XSec_dxdy_Ex();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double fx;
};

//.....................................................................................
//
// genie::utils::gsl::d2XSec_dxdy_Ey 
// A 1-D cross section function: d2xsec/dxdy = f(x)|(fixed:E,y) 
//
class d2XSec_dxdy_Ey: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dxdy_Ey(const XSecAlgorithmI * m, const Interaction * i, double x);
 ~d2XSec_dxdy_Ey();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fy;
};

//.....................................................................................
//
// genie::utils::gsl::d2XSec_dWdQ2_EW
// A 1-D cross section function: d2xsec/dWdQ2= f(Q2)|(fixed:E,W)
//
class d2XSec_dWdQ2_EW: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dWdQ2_EW( const XSecAlgorithmI * m, const Interaction * i, double W);
 ~d2XSec_dWdQ2_EW();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fW;
};

//.....................................................................................
//
// genie::utils::gsl::d2XSec_dWdQ2_EQ2
// A 1-D cross section function: d2xsec/dWdQ2= f(W)|(fixed:E,Q2)
//
class d2XSec_dWdQ2_EQ2: public ROOT::Math::IBaseFunctionOneDim
{
public:
  d2XSec_dWdQ2_EQ2(const XSecAlgorithmI * m, const Interaction * i, double Q2);
 ~d2XSec_dWdQ2_EQ2();

  // ROOT::Math::IBaseFunctionOneDim interface
  unsigned int                      NDim   (void)             const;
  double                            DoEval (double xin) const;
  ROOT::Math::IBaseFunctionOneDim * Clone  (void)             const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double                 fQ2;
};

//.....................................................................................

//.....................................................................................
//
// 
//
class d5XSecAR : public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d5XSecAR(const XSecAlgorithmI * m, const Interaction * i);
  ~d5XSecAR();
  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
  void SetFlip(bool b) { flip = b; }

private:
  const XSecAlgorithmI * fModel;
  const Interaction * fInteraction;
  bool flip;
};


//.....................................................................................
//
// genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi
// A 5-D cross section function (fixed E_nu)
//
class d5Xsec_dEldOmegaldOmegapi: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d5Xsec_dEldOmegaldOmegapi(const XSecAlgorithmI * m, const Interaction * i);
 ~d5Xsec_dEldOmegaldOmegapi();
        
  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
};

///.....................................................................................
///
/// genie::utils::gsl::d4Xsec_dEldThetaldOmegapi
/// A 4-D cross section function (fixed E_nu)
/// DANIEL - for the Alvarez-Russo cross-section
///
class d4Xsec_dEldThetaldOmegapi: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d4Xsec_dEldThetaldOmegapi(const XSecAlgorithmI * m, const Interaction * i);
 ~d4Xsec_dEldThetaldOmegapi();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
  
  double                              GetFactor()                 const;
  void                                SetFactor(double factor);

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  double fFactor;
};
///.....................................................................................
///
/// genie::utils::gsl::d3Xsec_dOmegaldThetapi
/// A 3-D cross section function (fixed E_nu)
/// Steve Dennis - for the Alvarez-Russo cross-section
///
class d3Xsec_dOmegaldThetapi: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3Xsec_dOmegaldThetapi(const XSecAlgorithmI * m, const Interaction * i);
 ~d3Xsec_dOmegaldThetapi();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  d3Xsec_dOmegaldThetapi            * Clone  (void)               const;
  
  // Specific to this class
  void SetE_lep (double E_lepton) const;
  // Yes, it's a const setter
  // Needed because DoEval must be const, but dXSec_dElep_AR::DoEval() must call this

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  mutable double fElep;
};
///.....................................................................................
///
/// genie::utils::gsl::dXSec_dElep_AR
/// A 1-D cross section function: dxsec/dElep
/// Used for Alvarez-Ruso coherent.
///
class dXSec_dElep_AR: public ROOT::Math::IBaseFunctionOneDim
{
public:
  dXSec_dElep_AR(const XSecAlgorithmI * m, const Interaction * i,
                 string gsl_nd_integrator_type, double gsl_relative_tolerance,
                 unsigned int max_n_calls);
  dXSec_dElep_AR() {};
 ~dXSec_dElep_AR();

  // ROOT::Math::IBaseFunctionOneDim interface
  dXSec_dElep_AR *                  Clone  (void)             const;
  double                            DoEval (double xin) const;

private:
  const XSecAlgorithmI * fModel;
  const Interaction *    fInteraction;
  
  const genie::utils::gsl::d3Xsec_dOmegaldThetapi * func;
  
  mutable ROOT::Math::IntegratorMultiDim integrator;
  
  double kine_min[3];
  double kine_max[3];
  
  string fGSLIntegratorType;
  double fGSLRelTol;
  unsigned int fGSLMaxCalls;
};
                  
///.....................................................................................
/// 
/// dXSec_Log_Wrapper
/// Redistributes variables over a range to a e^-x distribution.
/// Allows the integrator to use a logarithmic series of points while calling uniformly.
class dXSec_Log_Wrapper: public ROOT::Math::IBaseFunctionMultiDim
{
  public:
    dXSec_Log_Wrapper(const ROOT::Math::IBaseFunctionMultiDim * fn,
                      bool * ifLog, double * min, double * maxes);
   ~dXSec_Log_Wrapper();
  
    // ROOT::Math::IBaseFunctionMultiDim interface
    unsigned int                        NDim   (void)               const;
    double                              DoEval (const double * xin) const;
    ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
  
  private:
    const ROOT::Math::IBaseFunctionMultiDim * fFn;
    bool * fIfLog;
    double * fMins;
    double * fMaxes;
};
  
} // gsl   namespace
} // utils namespace
} // genie namespace

#endif   

