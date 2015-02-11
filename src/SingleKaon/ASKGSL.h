//_____________________________________________________________________________________
/*!

\brief      GSL wrappers for ASK model

\author     Chris Marshall and Martti Nirkko

\created    Feb 14, 2014

\cpright    Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//_____________________________________________________________________________________

#ifndef _ASK_GSL_WRAPPERS_H_
#define _ASK_GSL_WRAPPERS_H_

#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

namespace genie {

 class XSecAlgorithmI;
 class Interaction;

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

#endif   

