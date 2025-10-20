//____________________________________________________________________________
/*!

\class    genie::SmithMonizQELCCXSec

\brief    Computes the Quasi Elastic (QEL) cross section by Smith Moniz model. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research \n

          adapted from  fortran code provided by: \n

          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>
          Joint Institute for Nuclear Research \n

          Vadim Naumov <vnaumov@theor.jinr.ru>
          Joint Institute for Nuclear Research  \n

          based on code of: \n
          Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 05, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _SMITH_MONIZ_QEL_XSEC_H_
#define _SMITH_MONIZ_QEL_XSEC_H_

#include <Math/IFunction.h>

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"

namespace genie {



class SmithMonizQELCCXSec : public XSecIntegratorI {

public:
  SmithMonizQELCCXSec();
  SmithMonizQELCCXSec(string config);
  virtual ~SmithMonizQELCCXSec();

  //! XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  //! Overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

protected:
  string fGSLIntgType2D;  ///< name of GSL 2D numerical integrator
  double fGSLRelTol2D;    ///< required relative tolerance (error) for 2D integrator

private:
  mutable SmithMonizUtils * sm_utils;

  void LoadConfig (void);

};

//_____________________________________________________________________________________
//
// GSL wrappers
//
//_____________________________________________________________________________________

 namespace utils {
  namespace gsl   {

   class d2Xsec_dQ2dv: public ROOT::Math::IBaseFunctionMultiDim
   {
    public:
      d2Xsec_dQ2dv(const XSecAlgorithmI * m, const Interaction * i);
     ~d2Xsec_dQ2dv();
      // ROOT::Math::IBaseFunctionMultiDim interface
      unsigned int                        NDim   (void)               const;
      double                              DoEval (const double * xin) const;
      ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

    private:
      const XSecAlgorithmI * fModel;
      const Interaction *    fInteraction;
      mutable SmithMonizUtils * sm_utils;
      Range1D_t rQ2;
   };


  } // gsl   namespace
 } // utils namespace

}       // genie namespace
#endif  // _SMITH_MONIZ_QEL_XSEC_H_
