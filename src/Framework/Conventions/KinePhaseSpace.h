//____________________________________________________________________________
/*!

\class    genie::KinePhaseSpace

\brief    Enumeration of kinematical phase spaces

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 06, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _KINEMATIC_PHASE_SPACE_ENUM_H_
#define _KINEMATIC_PHASE_SPACE_ENUM_H_

#include <cassert>
#include <string>

using std::string;

namespace genie {
// Note: please attach new phase space enum element to the end of the list .
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
  kPSQD2fE,
  kPSlogQ2fE,
  kPSQ2fEW,
  kPSlogQ2fEW,
  kPSq2fE,
  kPSq2fEW,
  kPSWfE,
  kPSWfEQ2,
  kPSWfEq2,
  kPSWQ2fE,
  kPSWQD2fE,
  kPSW2Q2fE,
  kPSWlogQ2fE,
  kPSW2logQ2fE,
  kPSWq2fE,
  kPSW2q2fE,
  kPSxytfE,
  kPSQ2yfE,
  kPSlogQ2logyfE,
  kPSTlctl,
  kPSElOlOpifE,
  kPSElOlTpifE,
  kPSTkTlctl,
  kPSQ2vfE,
  kPSQELEvGen,// Phase space used by genie::QELEventGenerator for sampling kinematic variables
              // TODO: rename this value when the correct variables are identified
  kPSTAfE,
  kPSEgTlOgfE,
  kPSDMELEvGen, // Equivalent to kPSQELEvGen for Dark Matter scattering  
  kPSxQ2fE,
  kPSlog10xlog10Q2fE,
  kPSEDNufE, // Used for Dark Neutrinos, two body final state
  kPSn1n2fE,
  kPSn1n2n3fE,
  kPSWQ2ctpphipfE,
  kPSWQ2ctpfE,
  kPSQ2vpfE
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
      case(kPSxQ2fE)      : return "<{x,Q2}|E>";      break;
      case(kPSQ2fE)       : return "<{Q2}|E>";        break;
      case(kPSQD2fE)      : return "<{QD2}|E>";       break;
      case(kPSlogQ2fE)    : return "<{logQ2}|E>";     break;
      case(kPSQ2fEW)      : return "<{Q2}|E,W>";      break;
      case(kPSlogQ2fEW)   : return "<{logQ2}|E,W>";   break;
      case(kPSq2fE)       : return "<{q2}|E>";        break;
      case(kPSq2fEW)      : return "<{q2}|E,W>";      break;
      case(kPSWfE)        : return "<{W}|E>";         break;
      case(kPSWfEQ2)      : return "<{W}|E,Q2>";      break;
      case(kPSWfEq2)      : return "<{W}|E,q2>";      break;
      case(kPSWQ2fE)      : return "<{W,Q2}|E>";      break;
      case(kPSWQD2fE)     : return "<{W,QD2}|E>";     break;
      case(kPSW2Q2fE)     : return "<{W2,Q2}|E>";     break;
      case(kPSWlogQ2fE)   : return "<{W,logQ2}|E>";   break;
      case(kPSW2logQ2fE)  : return "<{W2,logQ2}|E>";  break;
      case(kPSWq2fE)      : return "<{W,q2}|E>";      break;
      case(kPSW2q2fE)     : return "<{W2,q2}|E>";     break;
      case(kPSxytfE)      : return "<{x,y,t}|E>";     break;
      case(kPSQ2yfE)      : return "<{Q2,y}|E>";      break;
      case(kPSlogQ2logyfE): return "<{Q2,y}|E>";      break;
      case(kPSTlctl)      : return "<{Tl,cos(theta_l)}|E>";     break;
      case(kPSElOlOpifE)  : return "<{El,Omega_l,Omega_pi}|E>"; break;
      case(kPSEgTlOgfE)   : return "<{Egamma,Theta_l,Omega_gamma}|E>"; break;
      case(kPSElOlTpifE)  : return "<{El,Omega_l,Theta_pi}|E>"; break;
      case(kPSTkTlctl)    : return "<{Tk,Tl,cos(theta_l)}|E>";  break;
      case(kPSQ2vfE)      : return "<{Q2,v}|E>"; break;
      case(kPSQ2vpfE)     : return "<{Q2,v,p}|E>"; break;
      // TODO: update this string when the appropriate kinematic variables are known
      case(kPSQELEvGen)   : return "<QELEvGen>"; break;
      case(kPSDMELEvGen)  : return "<DMELEvGen>"; break;
      case(kPSTAfE)       : return "<{TA}|E>";   break;
      case(kPSlog10xlog10Q2fE) : return "<{log10x,log10Q2}|E>"; break;
      case(kPSEDNufE)     : return "<{EDNu}|E>"; break;
      case(kPSn1n2fE)     : return "<{n1,n2}|E>"; break;
      case(kPSn1n2n3fE)   : return "<{n1,n2,n3}|E>"; break;
      case(kPSWQ2ctpphipfE): return "<{W, Q2, cost(theta_pion), phi_pion}|E>"; break;
      case(kPSWQ2ctpfE)    : return "<{W, Q2, cost(theta_pion)}|E>"; break;
    }
    return "** Undefined kinematic phase space **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINEMATIC_PHASE_SPACE_ENUM_H_
