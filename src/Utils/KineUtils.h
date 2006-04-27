//____________________________________________________________________________
/*!

\namespace  genie::utils::kinematics

\brief      Kinematical limits for DIS,QEL,RES

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    November 26, 2004

\cpright    Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
            All rights reserved.
            For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KINE_UTILS_H_
#define _KINE_UTILS_H_

#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "Interaction/Interaction.h"
#include "Utils/Range1.h"

namespace genie {
namespace utils {

namespace kinematics
{
  Range1D_t  KineRange             (const Interaction * const i, KineVar_t k);
  double     PhaseSpaceVolume      (const Interaction * const i, KinePhaseSpace_t ps);

  Range1D_t  WRange                (const Interaction * const i);
  Range1D_t  Q2Range               (const Interaction * const i);
  Range1D_t  q2Range               (const Interaction * const i);

  Range1D_t  Q2Range_W             (const Interaction * const i, Range1D_t rW);

  double     CalcQ2                (const Interaction * const i);
  double     CalcW                 (const Interaction * const i);

  void       ApplyCutsToKineLimits (Range1D_t & r, double min, double max);

  double     EnergyThreshold       (const Interaction * const i);

  bool       IsAboveCharmThreshold (const Interaction * const i, double mc);
  double     SlowRescalingVar      (const Interaction * const i, double mc);

} // kinematics namespace
} // utils namespace
} // genie namespace

#endif // _KINE_UTILS_H_
