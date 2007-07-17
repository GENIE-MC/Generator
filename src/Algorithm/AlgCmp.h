//____________________________________________________________________________
/*!

\class    genie::AlgCmp

\brief    Encapsulates an enumeration of possible algorithm comparisons

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 22, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _ALG_CMP_H_
#define _ALG_CMP_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EAlgCmp {

   kAlgCmpUnknown    = -1, 
   kAlgCmpIdentical, 
   kAlgCmpDiffConfig,
   kAlgCmpDiffAlg

} AlgCmp_t; 
  

class AlgCmp {

 public:

  static char * AsString(AlgCmp_t alg) {
     switch (alg) {
     case kAlgCmpIdentical:   return "Algorithm [same], configuration [same]";  break;
     case kAlgCmpDiffConfig:  return "Algorithm [same], configuration [diff]";  break;
     case kAlgCmpDiffAlg:     return "Algorithm [diff]";                        break;
     case kAlgCmpUnknown:     return "Undefined algorithm comparison result";   break;
     default:                 break;
     }
     return " ";
  }

};

}
#endif
