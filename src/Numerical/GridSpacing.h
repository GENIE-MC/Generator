//____________________________________________________________________________
/*!

\class    genie::GridSpacing

\brief    Encapsulates an enumeration of grid spacing methods

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 22, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _GRID_SPACING_H_
#define _GRID_SPACING_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EGridSpacing {

   kGSpUnknown    = -1,
   kGSpLinear,
   kGSpLoge,
   kGSpLog10,
   kGSpIrregular

} GridSpacing_t;


class GridSpacing {
 public:
  static char * AsString(GridSpacing_t gsp) {
     switch (gsp) {
     case kGSpLinear:    return "Grid axis spacing: uniform in x";        break;
     case kGSpLoge:      return "Grid axis spacing: uniform in ln(x)";    break;
     case kGSpLog10:     return "Grid axis spacing: uniform in log(x)";   break;
     case kGSpIrregular: return "Grid axis spacing: irregular";           break;
     case kGSpUnknown:   return "Undefined grid axis spacing method";     break;
     default:            break;
     }
     return " ";
  }
};

}       // namespace genie
#endif  // _GRID_SPACING_H_
