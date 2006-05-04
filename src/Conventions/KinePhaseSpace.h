//____________________________________________________________________________
/*!

\class    genie::KinePhaseSpace

\brief    Enumeration of kinematical phase spaces

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KINE_PHASE_SPACE_H_
#define _KINE_PHASE_SPACE_H_

#include <cassert>
#include <string>

using std::string;

namespace genie {

typedef enum EKinePhaseSpace {
  kPSNull = 0,
  kPSfE,
  kPSxfE,
  kPSlogxfE,
  kPSxfEy,
  kPSlogxfEy,
  kPSyfE,
  kPSlogyfE,
  kPSyfEx,
  kPSlogyfEx,
  kPSxyfE,
  kPSlogxlogyfE,
  kPSQ2fE,
  kPSlogQ2fE,
  kPSQ2fEW,
  kPSlogQ2fEW,
  kPSq2fE,
  kPSq2fEW,
  kPSWfE,
  kPSWfEQ2,
  kPSWfEq2,
  kPSWQ2fE,
  kPSWlogQ2fE,
  kPSWq2fE,
  kPSxytfE
} KinePhaseSpace_t;

class KinePhaseSpace
{
public:
  //__________________________________________________________________________
  static string AsString(KinePhaseSpace_t kps)
  {
    switch (kps) {

      case(kPSNull) : 
        return "** Undefined kinematic phase space **"; break;

      case(kPSfE)         : return "<|E>";            break;
      case(kPSxfE)        : return "<{x}|E>";         break;
      case(kPSlogxfE)     : return "<{logx}|E>";      break;
      case(kPSxfEy)       : return "<{x}|E,y>";       break;
      case(kPSlogxfEy)    : return "<{logx}|E,y>";    break;
      case(kPSyfE)        : return "<{y}|E>";         break;
      case(kPSlogyfE)     : return "<{logy}|E>";      break;
      case(kPSyfEx)       : return "<{y}|E,x>";       break;
      case(kPSlogyfEx)    : return "<{logy}|E,x>";    break;
      case(kPSlogxlogyfE) : return "<{logx,logy}|E>"; break;
      case(kPSxyfE)       : return "<{x,y}|E>";       break;
      case(kPSQ2fE)       : return "<{Q2}|E>";        break;
      case(kPSlogQ2fE)    : return "<{logQ2}|E>";     break;
      case(kPSQ2fEW)      : return "<{Q2}|E,W>";      break;
      case(kPSlogQ2fEW)   : return "<{logQ2}|E,W>";   break;
      case(kPSq2fE)       : return "<{q2}|E>";        break;
      case(kPSq2fEW)      : return "<{q2}|E,W>";      break;
      case(kPSWfE)        : return "<{W}|E>";         break;
      case(kPSWfEQ2)      : return "<{W}|E,Q2>";      break;
      case(kPSWfEq2)      : return "<{W}|E,q2>";      break;
      case(kPSWQ2fE)      : return "<{W,Q2}|E>";      break;
      case(kPSWlogQ2fE)   : return "<{W,logQ2}|E>";   break;
      case(kPSWq2fE)      : return "<{W,q2}|E>";      break;
      case(kPSxytfE)      : return "<{x,y,t}|E>";     break;
    }
    return "** Undefined kinematic phase space **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINE_PHASE_SPACE_H_
