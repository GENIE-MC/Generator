//____________________________________________________________________________
/*!

\class    genie::NievesQELCCXSec

\brief    Computes the Quasi Elastic (QEL) total cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Igor Kaorin JINR

\created  March 23, 2025

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NIEVES_QEL_XSEC_H_
#define _NIEVES_QEL_XSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"
#include "Physics/QuasiElastic/XSection/NievesQELCCPXSec.h"

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/AdaptiveIntegratorMultiDim.h>

namespace genie {   
class NievesQELCCXSec : public XSecIntegratorI {

public:

  NievesQELCCXSec(void);
  NievesQELCCXSec(std::string config);

  /// XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI* model, const Interaction* i) const;
  
  std::string Get1DimIntgType() const
  {
      return fGSL1DimIntgType;
  }
  
  unsigned int Get1DimMaxEval() const
  {
      return fGSL1DimMaxEval;
  }
  
  double Get1DimRelTol() const
  {
      return fGSL1DimRelTol;
  }

  /// Overload the Algorithm::Configure() methods to load private data
  /// members from configuration options
  void Configure(const Registry& config);
  void Configure(std::string config);

private:

  void LoadConfig (void);

  // XML configuration parameters
  std::string fGSLIntgType;
  double fGSLRelTol;
  unsigned int fGSLMaxEval;
  
  std::string fGSL1DimIntgType;
  double fGSL1DimRelTol;
  unsigned int fGSL1DimMaxEval;

};

class XSecAlgorithmI;
class Interaction;

namespace utils {
namespace gsl   {
//.....................................................................................
//
// A 3-D cross section function: d3xsec/dEldCosThetadR = f(El,CosTheta, R)|(fixed E)
//
class d3XSec_dElepdCosThetalepdR_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3XSec_dElepdCosThetalepdR_E(const XSecAlgorithmI * m, const Interaction * i, double Rmax);
 ~d3XSec_dElepdCosThetalepdR_E();

  // ROOT::Math::IBaseFunctionMultiDim interface
  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  const NievesQELCCPXSec * fXsec_model;
  const Interaction *    fInteraction;
  mutable SmithMonizUtils * sm_utils;
  Kinematics * fKinematics;
  double fRmax;
  double fEnu;
  double fml;
  double fml2;
};

} // gsl   namespace
} // utils namespace

} // genie namespace

#endif  // _NIEVES_QEL_XSEC_H_
