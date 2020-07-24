//____________________________________________________________________________
/*!

  \class      genie::utils::math::GTrace

  \brief      Simple math container for 4 x 4 complex objects 

  \author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
              University of Liverpool & STFC Rutherford Appleton Lab

              Marco Roda <mroda@liverpool.ac.uk>
              University of Liverpool

  \created    May 15, 2020

  \cpright    Copyright (c) 2003-2019, The GENIE Collaboration
              For the full text of the license visit http://copyright.genie-mc.org
              or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MATH_UTILS_GTRACE_H_
#define _MATH_UTILS_GTRACE_H_

#include "Framework/Numerical/MathUtils.h"

namespace genie {
  namespace utils {

    namespace math  {

      class GTrace : public ComplexMatrix<4> {

      public:

        const GTrace &operator+=( const GTrace & ) ;
        const GTrace &operator*=( std::complex<double> ) ;
        const GTrace &operator*=( double c ) { return (*this) *= std::complex<double>(c, 0) ; } ;

        GTrace Conj( const GTrace & M ) ;

        //GTrace operator*(GTrace &);
        //GTrace operator+(double);
        //GTrace operator*(double);
	
      };

      class TraceContraction : public std::pair<GTrace, GTrace> {
      public:

    	  TraceContraction( const GTrace & a, const GTrace & b ) :
    	  std::pair<GTrace, GTrace>(a, b) { ; }

    	  std::complex<double> operator()(unsigned short i, unsigned short j,
    			                          unsigned short m, unsigned short n) {
          return first[i][j] * std::conj(second[m][n]) ; }
      };

    } // math  namespace
  } // utils namespace
} // genie namespace

#endif // _MATH_UTILS_GTRACE_H_ 
