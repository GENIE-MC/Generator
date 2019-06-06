//____________________________________________________________________________
/*!

\class    genie::AlamSimoAtharVacasSKXSec

\brief    A cross-section integrator and GSL interface for the 
          M. Rafi Alam, I. Ruiz Simo, M. Sajjad Athar and M.J. Vicente Vacas
          single-Kaon production model.
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Chris Marshall and Martti Nirkko

\created  March 20, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_XSEC_H_
#define _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_XSEC_H_

#include <Math/Integrator.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

namespace genie {

class XSecAlgorithmI;
class Interaction;

class AlamSimoAtharVacasSKXSec : public XSecIntegratorI {
public:
  AlamSimoAtharVacasSKXSec();
  AlamSimoAtharVacasSKXSec(string config);
  virtual ~AlamSimoAtharVacasSKXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

protected:
  bool fSplitIntegral;

private:
  void LoadConfig (void);
};

//_____________________________________________________________________________________
// 
// GSL wrappers 
//
//_____________________________________________________________________________________

 namespace utils {
  namespace gsl   {

   class d3Xsec_dTldTkdCosThetal: public ROOT::Math::IBaseFunctionMultiDim
   {
    public:
      d3Xsec_dTldTkdCosThetal(const XSecAlgorithmI * m, const Interaction * i);
     ~d3Xsec_dTldTkdCosThetal();
      // ROOT::Math::IBaseFunctionMultiDim interface
      unsigned int                        NDim   (void)               const;
      double                              DoEval (const double * xin) const;
      ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;
    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
   };

  } // gsl   namespace
 } // utils namespace

} // genie namespace

#endif  // _ALAM_SIMO_ATHAR_VACAS_SINGLE_KAON_XSEC_H_
