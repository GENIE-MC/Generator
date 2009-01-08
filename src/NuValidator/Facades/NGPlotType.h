//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGPlotType

\brief    NeuGEN's Plot type enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NG_PLOT_TYPE_H_
#define _NG_PLOT_TYPE_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGPlotType {

  e_XSec = 0,
  e_DiffXSec,
  e_SF,
  e_UndefinedPlotType

} NGPlotType_t;

class NGPlotType {

  public:

     virtual ~NGPlotType() { }

     static const char * AsString(NGPlotType_t plot_type) 
     {  
       switch(plot_type) {
         case e_XSec:       return "cross section";              break;
         case e_DiffXSec:   return "differential cross section"; break;
         case e_SF:         return "structure function";         break;
  
         case e_UndefinedPlotType:
         default:            
            return "Unknown plot type"; break;
       }
       return "Unknown plot type"; 
     }
     
     static NGPlotType_t GetPlotTypeFromCode(int code)
     {
        if      (code ==  0) return e_XSec;
        else if (code ==  1) return e_DiffXSec;
        else if (code ==  2) return e_SF;
        else                 return e_UndefinedPlotType;
     }

ClassDef(NGPlotType, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace


#endif

