//____________________________________________________________________________
/*!

\class    genie::AlgStatus

\brief    Encapsulates an enumeration of possible algorithm execution states.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004
 
*/
//____________________________________________________________________________

#ifndef _ALG_STATUS_H_
#define _ALG_STATUS_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EAlgStatus {

   kAlgUndefinedStatus    = -1, 
   kAlgFail,
   kAlgSuccess,
   kAlgNotConverged

} AlgStatus_t; 
  

class AlgStatus {

 public:

  static char * AsString(AlgStatus_t alg) {
     switch (alg) {
     case kAlgFail:                 return "FAIL";             break;
     case kAlgSuccess:              return "SUCCESS";          break;
     case kAlgUndefinedStatus:      return "UNDEFINED STATUS"; break;
     case kAlgNotConverged:         return "NOT CONVERGED";    break;
     default:                       break;
     }
     return " ";
  }

};

}
#endif
