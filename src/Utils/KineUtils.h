//____________________________________________________________________________
/*!

\namespace  genie::utils::kinematics

\brief      Kinematical utilities

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
  //! methods used for computing phase space volumes
  double PhaseSpaceVolume (const Interaction * const i, KinePhaseSpace_t ps);

  //! methods used for computing the Jacobians associated with phase space transformations
  double Jacobian         (const Interaction * const i, KinePhaseSpace_t f, KinePhaseSpace_t t);
  bool   TransformMatched (KinePhaseSpace_t ia, KinePhaseSpace_t ib,
                               KinePhaseSpace_t a, KinePhaseSpace_t b, bool & fwd);

  //! methods used for figuring out the physical range of kinematical variables
  double     EnergyThreshold (const Interaction * const i);
  Range1D_t  KineRange       (const Interaction * const i, KineVar_t k);
  Range1D_t  WRange          (const Interaction * const i);
  Range1D_t  Q2Range         (const Interaction * const i);
  Range1D_t  q2Range         (const Interaction * const i);
  Range1D_t  Q2Range_W       (const Interaction * const i, Range1D_t rW);
  void       MinXY           (const Interaction * const i, double & x, double & y);

  //! methods used to apply cuts to kinematical limits
  void ApplyCutsToKineLimits (Range1D_t & r, double min, double max);

  //! kinematical variable transforms
  double QD2toQ2 (double QD2);
  double Q2toQD2 (double Q2);
  void   WQ2toXY (double Ev, double M, double W, double Q2, double & x, double & y);
  void   XYtoWQ2 (double Ev, double M, double & W, double & Q2, double x, double y);
  double XYtoW   (double Ev, double M, double x, double y);
  double XYtoQ2  (double Ev, double M, double x, double y);

  double CalcQ2  (const Interaction * const i);
  double CalcW   (const Interaction * const i);

  //! charm threshold and slow rescaling variable
  bool   IsAboveCharmThreshold (const Interaction * const i, double mc);
  double SlowRescalingVar      (const Interaction * const i, double mc);

  //! Functions used to define differential cross section distribution envelopes 
  //! for importance sampling in kinematical selection modules
  double RESImportanceSamplingEnvelope(double * x, double * par);
  double DISImportanceSamplingEnvelope(double * x, double * par);

} // kinematics namespace
} // utils namespace
} // genie namespace

#endif // _KINE_UTILS_H_
