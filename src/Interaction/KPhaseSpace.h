//____________________________________________________________________________
/*!

\class    genie::KPhaseSpace

\brief    Kinematical phase space 

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _KINEMATIC_PHASE_SPACE_H_
#define _KINEMATIC_PHASE_SPACE_H_

#include <cassert>

#include <TObject.h>

#include "Conventions/KineVar.h"
//#include "Interaction/KPhaseSpaceCut.h"
#include "Utils/Range1.h"

namespace genie {

class Interaction;

class KPhaseSpace : public TObject {

public:
  KPhaseSpace (void);
  KPhaseSpace (const Interaction * in);
 ~KPhaseSpace (void);

  void UseInteraction(const Interaction * in);

  //-- Energy threshold
  double Threshold(void) const;

  //-- Checks whether the interaction is above the energy threshold
  bool IsAboveThreshold(void) const;

  //-- Check whether the running kinematics for the specified interaction
  //-- are allowed. Takes into account all applied cuts.
  bool IsAllowed (void) const;

  //-- Return the kinematical variable limits
  Range1D_t Limits  (KineVar_t kvar) const;
  double    Minimum (KineVar_t kvar) const;
  double    Maximum (KineVar_t kvar) const;

  //-- Apply kinematical cuts
  //void ApplyCut (const KPhaseSpaceCut & cut);
  //void ApplyCut (KineVar_t kvar, Range1D_t narrower_range);

  Range1D_t  WLim    (void) const;
  Range1D_t  Q2Lim_W (void) const;
  Range1D_t  q2Lim_W (void) const;
  Range1D_t  Q2Lim   (void) const;
  Range1D_t  q2Lim   (void) const;
  Range1D_t  XLim    (void) const;
  Range1D_t  YLim    (void) const;
  Range1D_t  YLim_X  (void) const;

private:
  void Init(void);
  const Interaction * fInteraction;

ClassDef(KPhaseSpace,1)
};

}      // genie namespace
#endif // _KINE_PHASE_SPACE_H_

