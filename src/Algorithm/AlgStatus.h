//____________________________________________________________________________
/*!

\class    genie::AlgStatus

\brief    Encapsulates an enumeration of possible algorithm execution states.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
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
   kAlgSuccess

} AlgStatus_t; 
  

class AlgStatus {

public:

  static char * AsString(AlgStatus_t alg) {
     switch (alg) {
     case kAlgFail:               return "Algorithm failed";           break;
     case kAlgSuccess:            return "Algorithm run successfully"; break;
     case kAlgUndefinedStatus:    return "Undefined alg status";       break;
     default:                     break;
     }
     return " ";
  }

};

}
#endif
