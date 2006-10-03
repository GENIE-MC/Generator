//____________________________________________________________________________
/*!

\class    genie::IntranukeMode

\brief    An enumeration of intranuke modes

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.

\created  October 3, 2006

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE. 
*/
//____________________________________________________________________________

#ifndef _INTRANUKE_MODE_H_
#define _INTRANUKE_MODE_H_

#include <string>

using std::string;

namespace genie {

typedef enum EINukeMode {

   kIMdUndefined = -1, 
   kIMdhN,
   kIMdhA

} INukeMode_t;   

class INukeMode {

public:
  //__________________________________________________________________________
  static string AsString(INukeMode_t mode) {
     switch (mode) {
     case kIMdUndefined:  return "** Undefined Intranuke mode **"; break;
     case kIMdhN:         return "hadron+nucleon (h+N)";  break;
     case kIMdhA:         return "hadron+nucleus (h+A)";  break;
     default:             break;
     }
     return "** Undefined Intranuke mode ** ";
  }
  //__________________________________________________________________________
};

}      // genie
#endif // _INTRANUKE_MODE_H_
