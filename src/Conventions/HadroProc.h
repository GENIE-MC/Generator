//____________________________________________________________________________
/*!

\class    genie::HadroProc

\brief    An enumeration of possible hadron interaction processes taken into 
          account by intranuclear rescattering event generator modules

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 1, 2005
 
*/
//____________________________________________________________________________

#ifndef _HADRONIC_PROCESSES_H_
#define _HADRONIC_PROCESSES_H_

#include <string>

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

using std::string;

namespace genie {

typedef enum EHadroProc {

   kHPcUndefined      = -1, 
   kHPcChargeExchange,
   kHPcInelastic,
   kHPcAbsorption,
   kHPcElastic

} HadroProc_t;   

class HadroProc {

public:
  //__________________________________________________________________________
  static string AsString(HadroProc_t p) {
     switch (p) {
     case kHPcUndefined:       return "** Undefined hadronic process **";break;
     case kHPcChargeExchange:  return "Charge Exchange";  break;
     case kHPcInelastic:       return "Inelastic";        break;
     case kHPcAbsorption:      return "Absorption";       break;
     case kHPcElastic:         return "Elastic";          break;
     default:                   break;
     }
     return "** Undefined hadronic process ** ";
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _HADRONIC_PROCESSES_H_
