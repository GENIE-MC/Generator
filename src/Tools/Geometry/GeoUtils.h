//____________________________________________________________________________
/*!

\namespace  genie::utils::geometry

\brief      Geometry utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    March 26, 2009

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GEO_UTILS_H_
#define _GEO_UTILS_H_

#include <string>

class TGeoVolume;

using std::string;

namespace genie {
namespace utils {

namespace geometry   
{
  void RecursiveExhaust(TGeoVolume *topvol, string volnames, bool exhaust);

} // geometry namespace
} // utils namespace
} // genie namespace

#endif // _GEO_UTILS_H_
