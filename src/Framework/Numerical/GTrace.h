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

        GTrace operator+(GTrace &);
        GTrace operator*(GTrace &);
	
      };

    } // math  namespace
  } // utils namespace
} // genie namespace

#endif // _MATH_UTILS_GTRACE_H_ 
