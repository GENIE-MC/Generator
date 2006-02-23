//____________________________________________________________________________
/*!

\namespace  genie::utils::kinematics

\brief      Kinematical limits for DIS,QEL,RES

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    Novemner 26, 2004

*/
//____________________________________________________________________________

#ifndef _KINE_UTILS_H_
#define _KINE_UTILS_H_

#include "Interaction/Interaction.h"
#include "Utils/Range1.h"

namespace genie {
namespace utils {

namespace kinematics
{
  Range1D_t  WRange      (const Interaction * const interaction);

  Range1D_t  Q2Range     (const Interaction * const interaction);
  Range1D_t  q2Range     (const Interaction * const interaction);

  Range1D_t  Q2Range_W   (const Interaction * const interaction);
  Range1D_t  Q2Range_W   (const Interaction * const interaction, Range1D_t rW);

  Range1D_t  Q2Range_xy  (const Interaction * const interaction);
  Range1D_t  Q2Range_M   (const Interaction * const interaction);

  Range1D_t  q2Range_W   (const Interaction * const interaction);
  Range1D_t  q2Range_xy  (const Interaction * const interaction);
  Range1D_t  q2Range_M   (const Interaction * const interaction);

  double     EnergyThreshold       (const Interaction * const interaction);
  bool       IsAboveCharmThreshold (const Interaction * const interaction, double mc);
  void       ApplyCutsToKineLimits (Range1D_t & r, double min, double max);
  double     CalcQ2                (const Interaction * const interaction);
  double     SlowRescalingVar      (const Interaction * const interaction, double mc);

} // kinematics namespace
} // utils namespace
} // genie namespace

#endif // _KINE_UTILS_H_
