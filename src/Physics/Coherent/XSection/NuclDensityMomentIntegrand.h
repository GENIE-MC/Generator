//____________________________________________________________________________
/*!
\class    genie::utils::gsl::wrap::NuclDensityMomentIntegrand

\brief    Integrand for the calculation of the k^th nuclear density moment: 
	  \int_{0}^{\infinity} \rho_{A}(r) r^k d^{3}r
          where \rho_{A}(r) is the nuclear density for a nucleus with atomic
          mass number A.
         
\author   Costas Andreopoulos <costas.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory      

\created  July 12, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NUCL_DENSITY_MOMENT_INTEGRAND_H_
#define _NUCL_DENSITY_MOMENT_INTEGRAND_H_

#include <Math/IFunction.h>

namespace genie {
 namespace utils {
  namespace gsl   {
   namespace wrap   {

    class NuclDensityMomentIntegrand : public ROOT::Math::IBaseFunctionOneDim
    {
     public:
       NuclDensityMomentIntegrand(int A, int k /*, would be good to have options on \rho */);
      ~NuclDensityMomentIntegrand();
       // ROOT::Math::IBaseFunctionOneDim interface
       unsigned int                      NDim   (void)       const;
       double                            DoEval (double xin) const;
       ROOT::Math::IBaseFunctionOneDim * Clone  (void)       const;
     private:
       int fA; // atomic mass
       int fK; // which nuclear density moment
    };

   } // wrap namespace
  } // gsl namespace
 } // utils namespace
} // genie namespace


#endif  // _NUCL_DENSITY_MOMENT_INTEGRAND_H_
