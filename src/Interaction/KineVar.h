//____________________________________________________________________________
/*!

\class    genie::KineVar

\brief    Enumeration of kinematic variables

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

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
  kKVW

} KineVar_t;

class KineVar
{
public:

  //__________________________________________________________________________
  static string AsString(KineVar_t kv)
  {
    switch (kv) {

      case(kKVNull) : return "[Undefined kinematic variable]"; break;
      case(kKVx)    : return "[Bjorken x]";                    break;
      case(kKVy)    : return "[Inelasticity y]";               break;
      case(kKVQ2)   : return "[Momentum transfer Q2 (>0)]";    break;
      case(kKVq2)   : return "[Momentum transfer q2 (<0)]";    break;
      case(kKVW)    : return "[Hadronic invariant mass W]";    break;
      default       : return "[Unknown kinematic variable]";   break;
    }
    return "[Unknown kinematic variable]";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINE_VAR_H_
