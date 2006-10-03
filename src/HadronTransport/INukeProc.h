//____________________________________________________________________________
/*!

\class    genie::INukeProc

\brief    An enumeration of possible hadron interaction processes taken into 
          account by the INTRANUKE intranuclear rescattering module.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.

\created  November 1, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_PROCESSES_H_
#define _INTRANUKE_PROCESSES_H_

#include <string>

using std::string;

namespace genie {

typedef enum EINukeProc {

   kIPcUndefined      = -1, 
   kIPcChargeExchange,
   kIPcInelastic,
   kIPcAbsorption,
   kIPcElastic

} INukeProc_t;   

class INukeProc {

public:
  //__________________________________________________________________________
  static string AsString(INukeProc_t proc) {
     switch (proc) {
     case kIPcUndefined:       return "** Undefined Intranuke process **";break;
     case kIPcChargeExchange:  return "Charge Exchange";  break;
     case kIPcInelastic:       return "Inelastic";        break;
     case kIPcAbsorption:      return "Absorption";       break;
     case kIPcElastic:         return "Elastic";          break;
     default:                  break;
     }
     return "** Undefined Intranuke process ** ";
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_PROCESSES_H_
