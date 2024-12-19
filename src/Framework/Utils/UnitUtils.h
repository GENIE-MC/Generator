//____________________________________________________________________________
/*!

\namespace  genie::utils::units

\brief      Simple unit system utilities

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    May 06, 2004

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _UNIT_UTILS_H_
#define _UNIT_UTILS_H_

#include <string>
using std::string;

namespace genie {
namespace utils {

namespace units {

  double UnitFromString(string u);

} // namespace units
} // namespace utils
} // namespace genie

#endif // _UNIT_UTILS_H_
