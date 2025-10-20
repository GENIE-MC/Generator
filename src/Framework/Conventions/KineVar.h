//____________________________________________________________________________
/*!

\class    genie::KineVar

\brief    Enumeration of kinematic variables

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 06, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org      
*/
//____________________________________________________________________________

#ifndef _KINEMATIC_VAR_ENUM_H_
#define _KINEMATIC_VAR_ENUM_H_

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
  kKVTk,
  kKVTl,
  kKVctl,
  kKVphikq,
  kKVSelx,
  kKVSely,
  kKVSelQ2,
  kKVSelq2,
  kKVSelW,
  kKVSelt,
  kKVSelTk,
  kKVSelTl,
  kKVSelctl,
  kKVSelphikq,
  kKVSelRad,
  kKVPn,
  kKVv,
  kKVphip,
  kKVctp,
  kKVSelv,
  kKVQ0,
  kKVQ3,
  kKVSelQ0,
  kKVSelQ3,
  kKVn1,
  kKVn2,
  kKVn3,
  kKVSelphip,
  kKVSelctp,
  // put all new enum names right before this line
  // do not change any previous ordering (neither insert nor delete)
  kNumOfKineVar

} KineVar_t;

class KineVar
{
public:

  //__________________________________________________________________________
  static string AsString(KineVar_t kv)
  {
    switch (kv) {
      case(kKVNull)    : return "** Undefined kinematic variable **";    break;
      case(kKVx)       : return " *Running* Bjorken x";                  break;
      case(kKVy)       : return " *Running* Inelasticity y";             break;
      case(kKVQ2)      : return " *Running* Momentum transfer Q2 (>0)";  break;
      case(kKVq2)      : return " *Running* Momentum transfer q2 (<0)";  break;
      case(kKVW)       : return " *Running* Hadronic invariant mass W";  break;
      case(kKVt)       : return " *Running* COH 4p transfer to nucleus"; break;
      case(kKVTk)      : return " *Running* meson kinetic energy";       break;
      case(kKVTl)      : return " *Running* lepton kinetic energy";      break;
      case(kKVctl)     : return " *Running* cosine of lepton theta";     break;
      case(kKVphikq)   : return " *Running* ASK phi kq";                 break;
      case(kKVSelx)    : return "*Selected* Bjorken x";                  break;
      case(kKVSely)    : return "*Selected* Inelasticity y";             break;
      case(kKVSelQ2)   : return "*Selected* Momentum transfer Q2 (>0)";  break;
      case(kKVSelq2)   : return "*Selected* Momentum transfer q2 (<0)";  break;
      case(kKVSelW)    : return "*Selected* Hadronic invariant mass W";  break;
      case(kKVSelt)    : return "*Selected* COH 4p transfer to nucleus"; break;
      case(kKVSelTk)   : return "*Selected* ASK kaon kinetic energy";    break;
      case(kKVSelTl)   : return "*Selected* ASK lepton kinetic energy";  break;
      case(kKVSelctl)  : return "*Selected* ASK cosine lepton theta";    break;
      case(kKVSelphikq): return "*Selected* ASK phi kq";                 break;
      case(kKVSelRad)  : return "*Selected* Struck particle position";   break;
      case(kKVPn)      : return " *Running* Hit nucleon momentum";       break;
      case(kKVv)       : return " *Running* Energy transfer";            break;
      case(kKVphip)    : return " *Running* SPP final pion phi";          break;
      case(kKVctp)     : return " *Running* SPP cosine final pion theta"; break;
      case(kKVSelv)    : return "*Selected* Energy transfer";            break;
      case(kKVQ0)      : return " *Running* Energy transfer (Q0)" ;      break;
      case(kKVQ3)      : return " *Running* three momentum transfer" ;   break;
      case(kKVSelQ0)   : return "*Selected* energy transfer (Q0) " ;     break;
      case(kKVSelQ3)   : return "*Selected* three momentum transfer" ;   break;
      case(kKVn1)      : return " *Running* Normalized variable n1"  ;   break;
      case(kKVn2)      : return " *Running* Normalized variable n2"  ;   break;
      case(kKVn3)      : return " *Running* Normalized variable n3"  ;   break;
      case(kKVSelphip) : return "*Selected* SPP final pion phi";          break;
      case(kKVSelctp)  : return "*Selected* SPP cosine final pion theta"; break;
      default          : return "** Unknown kinematic variable **";      break;
    }
    return "** Unknown kinematic variable **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINEMATIC_VAR_ENUM_H_
