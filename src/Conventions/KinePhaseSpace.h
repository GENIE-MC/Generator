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
        return "** Undefined kinematic phase space **"; 
        break;
      case(kPSfE) : 
        return "0-D Phase Space: |E"; 
        break;
      case(kPSxfE) : 
        return "1-D Phase Space: {x}|E"; 
        break;
      case(kPSlogxfE) : 
        return "1-D Phase Space: {logx}|E"; 
        break;
      case(kPSxfEy) : 
        return "1-D Phase Space: {x}|E,y"; 
        break;
      case(kPSlogxfEy) : 
        return "1-D Phase Space: {logx}|E,y";  
        break;
      case(kPSyfE) : 
        return "1-D Phase Space: {y}|E"; 
        break;
      case(kPSlogyfE) : 
        return "1-D Phase Space: {logy}|E";                
        break;
      case(kPSyfEx) : 
        return "1-D Phase Space: {y}|E,x";              
        break;
      case(kPSlogyfEx) : 
        return "1-D Phase Space: {logy}|E,x";  
        break;
      case(kPSxyfE) : 
        return "2-D Phase Space: {x,y}|E";  
        break;
      case(kPSQ2fE) :
        return "1-D Phase Space: {Q2}|E";  
        break;
      case(kPSlogQ2fE) :
        return "1-D Phase Space: {logQ2}|E";  
        break;
      case(kPSQ2fEW) : 
        return "1-D Phase Space: {Q2}|E,W";  
        break;
      case(kPSlogQ2fEW) : 
        return "1-D Phase Space: {logQ2}|E,W";  
        break;
      case(kPSq2fE) :
        return "1-D Phase Space: {q2}|E";   
        break;
      case(kPSq2fEW) : 
        return "1-D Phase Space: {q2}|E,W";  
        break;
      case(kPSWfE) :
        return "1-D Phase Space: {W}|E";  
        break;
      case(kPSWfEQ2) : 
        return "1-D Phase Space: {W}|E,Q2";
        break;
      case(kPSWfEq2) : 
        return "1-D Phase Space: {W}|E,q2";
        break;
      case(kPSWQ2fE) :
        return "2-D Phase Space: {W,Q2}|E"; 
        break;
      case(kPSWq2fE) :
        return "2-D Phase Space: {W,q2}|E";  
        break;
      case(kPSxytfE) : 
        return "3-D Phase Space: {x,y,t}|E";  
        break;
    }
    return "** Undefined kinematic phase space **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINE_PHASE_SPACE_H_
