//____________________________________________________________________________
/*!

\class    genie::InukeInt

\brief    Encapsulates an enumeration of possible hadron re-interactions 
          taken into account by Intranuke EventRecordVisitorI

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 1, 2005
 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_INTERACTION_H_
#define _INTRANUKE_INTERACTION_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EInukeInt {

   kInukiUndefined      = -1, 
   kInukiChargeExchange,
   kInukiInelastic,
   kInukiAbsorption,
   kInukiElastic

} InukeInt_t; 
  

class InukeInt {

 public:

  static char * AsString(InukeInt_t inuci) {
     switch (inuci) {
     case kInukiUndefined:       return "Undefined intranuclear interaction";  break;
     case kInukiChargeExchange:  return "Charge Exchange";  break;
     case kInukiInelastic:       return "Inelastic";        break;
     case kInukiAbsorption:      return "Absorption";       break;
     case kInukiElastic:         return "Elastic";          break;
     default:                    break;
     }
     return " ";
  }

};

}
#endif
