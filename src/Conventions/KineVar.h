//____________________________________________________________________________
/*!

\class    genie::KineVar

\brief    Enumeration of kinematic variables

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KINE_VAR_H_
#define _KINE_VAR_H_

#include <cassert>
#include <string>

using std::string;

namespace genie {

typedef enum EKineVar {

  kKVNull = 0,
  kKVx,
  kKVy,
  kKVQ2,
  kKVq2,
  kKVW,
  kKVt,
  kKVSelx,
  kKVSely,
  kKVSelQ2,
  kKVSelq2,
  kKVSelW,
  kKVSelt

} KineVar_t;

class KineVar
{
public:

  //__________________________________________________________________________
  static string AsString(KineVar_t kv)
  {
    switch (kv) {
      case(kKVNull) : return "** Undefined kinematic variable **";    break;
      case(kKVx)    : return " *Running* Bjorken x";                  break;
      case(kKVy)    : return " *Running* Inelasticity y";             break;
      case(kKVQ2)   : return " *Running* Momentum transfer Q2 (>0)";  break;
      case(kKVq2)   : return " *Running* Momentum transfer q2 (<0)";  break;
      case(kKVW)    : return " *Running* Hadronic invariant mass W";  break;
      case(kKVt)    : return " *Running* COH 4p transfer to nucleus"; break;
      case(kKVSelx) : return "*Selected* Bjorken x";                  break;
      case(kKVSely) : return "*Selected* Inelasticity y";             break;
      case(kKVSelQ2): return "*Selected* Momentum transfer Q2 (>0)";  break;
      case(kKVSelq2): return "*Selected* Momentum transfer q2 (<0)";  break;
      case(kKVSelW) : return "*Selected* Hadronic invariant mass W";  break;
      case(kKVSelt) : return "*Selected* COH 4p transfer to nucleus"; break;
      default       : return "** Unknown kinematic variable **";      break;
    }
    return "** Unknown kinematic variable **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINE_VAR_H_
