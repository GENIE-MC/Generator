//____________________________________________________________________________
/*!

\class    genie::IntranukeMode

\brief    An enumeration of intranuke modes

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.

\created  October 3, 2006

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_MODE_H_
#define _INTRANUKE_MODE_H_

#include <string>

using std::string;

namespace genie {

typedef enum EINukeMode {

   kIMdUndefined = -1, 
   kIMdHN,
   kIMdHA

} INukeMode_t;   

class INukeMode {

public:
  //__________________________________________________________________________
  static string AsString(INukeMode_t mode) {
     switch (mode) {
     case kIMdUndefined:  return "** Undefined Intranuke mode **"; break;
     case kIMdHN:         return "hN Intranuke mode";  break;
     case kIMdHA:         return "hA Intranuke mode";  break;
     default:             break;
     }
     return "** Undefined Intranuke mode ** ";
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_MODE_H_
