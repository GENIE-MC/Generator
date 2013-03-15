//____________________________________________________________________________
/*!

\class    genie::GridSpacing

\brief    Encapsulates an enumeration of grid spacing methods

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 22, 2005

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
  static const char * AsString(GridSpacing_t gsp) {
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
