//____________________________________________________________________________
/*!

\class    genie::IntranukeProc

\brief    Encapsulates an enumeration of possible hadron re-interaction
          processes taken into account by Intranuke EventRecordVisitorI

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 1, 2005
 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_PROCESSES_H_
#define _INTRANUKE_PROCESSES_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum EINukeProc {

   kINukUndefined      = -1, 
   kINukChargeExchange,
   kINukInelastic,
   kINukAbsorption,
   kINukElastic

} INukeProc_t;   

class INukeProc {

 public:

  static char * AsString(INukeProc_t p) {
     switch (p) {
     case kINukUndefined:       return "Undefined intranuclear interaction";  break;
     case kINukChargeExchange:  return "Charge Exchange";  break;
     case kINukInelastic:       return "Inelastic";        break;
     case kINukAbsorption:      return "Absorption";       break;
     case kINukElastic:         return "Elastic";          break;
     default:                   break;
     }
     return " ";
  }

};

}
#endif
