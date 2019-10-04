//____________________________________________________________________________
/*!

\namespace  genie::utils::gsl

\brief      Simple utilities for integrating GSL in the GENIE framework

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 06, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GSL_UTILS_H_
#define _GSL_UTILS_H_

#include <Math/AllIntegrationTypes.h>

namespace genie {
namespace utils {
namespace gsl   {

  ROOT::Math::IntegrationOneDim::Type   
       Integration1DimTypeFromString (string type);

  ROOT::Math::IntegrationMultiDim::Type 
       IntegrationNDimTypeFromString (string type);

} // namespace gsl
} // namespace utils
} // namespace genie

#endif // _GSL_UTILS_H_

