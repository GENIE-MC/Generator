//____________________________________________________________________________
/*!

\namespace  genie::utils::gsl

\brief      Simple utilities for integrating GSL in the GENIE framework

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    May 06, 2004

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org            
*/
//____________________________________________________________________________

#ifndef _GSL_UTILS_H_
#define _GSL_UTILS_H_

#include <Math/AllIntegrationTypes.h>
#include <string>
using std::string;

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
