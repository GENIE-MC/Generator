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
  kPSx,
  kPSy,
  kPSxy,
  kPSQ2,
  kPSq2,
  kPSW,
  kPSWQ2,
  kPSWq2,
  kPSxyt
} KinePhaseSpace_t;

class KinePhaseSpace
{
public:
  //__________________________________________________________________________
  static string AsString(KinePhaseSpace_t kps)
  {
    switch (kps) {
      case(kPSNull) : return "** Undefined kinematic phase space **"; break;
      case(kPSx)    : return "1-D Phase Space: {x}";                  break;
      case(kPSy)    : return "1-D Phase Space: {y}";                  break;
      case(kPSxy)   : return "2-D Phase Space: {x,y}";                break;
      case(kPSQ2)   : return "1-D Phase Space: {Q2}";                 break;
      case(kPSq2)   : return "1-D Phase Space: {q2}";                 break;
      case(kPSW)    : return "1-D Phase Space: {W}";                  break;
      case(kPSWQ2)  : return "2-D Phase Space: {W,Q2}";               break;
      case(kPSWq2)  : return "2-D Phase Space: {W,q2}";               break;
      case(kPSxyt)  : return "3-D Phase Space: {x,y,t}";              break;
    }
    return "** Undefined kinematic phase space **";
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _KINE_PHASE_SPACE_H_
