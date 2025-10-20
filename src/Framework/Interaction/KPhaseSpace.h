//____________________________________________________________________________
/*!

\class    genie::KPhaseSpace

\brief    Kinematical phase space

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 06, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _KINEMATIC_PHASE_SPACE_H_
#define _KINEMATIC_PHASE_SPACE_H_

#include <cassert>

#include <TObject.h>

#include "Framework/Conventions/KineVar.h"
//#include "Interaction/KPhaseSpaceCut.h"
#include "Framework/Utils/Range1.h"

namespace genie {

class Interaction;

class KPhaseSpace : public TObject {

public:
  KPhaseSpace (void);
  KPhaseSpace (const Interaction * in);
 ~KPhaseSpace (void);

  void UseInteraction(const Interaction * in);

  //! Energy threshold
  double Threshold(void) const;
  double Threshold_SPP_iso(void) const;   ///< Energy limit for resonance single pion production on isoscalar nucleon

  //! Checks whether the interaction is above the energy threshold
  bool IsAboveThreshold(void) const;

  //! Check whether the current kinematics is in the allowed phase space
  bool IsAllowed (void) const;

  //! Return the kinematical variable limits
  Range1D_t  Limits  (KineVar_t kvar) const;
  double     Minimum (KineVar_t kvar) const;
  double     Maximum (KineVar_t kvar) const;

  Range1D_t  WLim    (void) const;  ///< W  limits
  Range1D_t  Q2Lim_W (void) const;  ///< Q2 limits @ fixed W
  Range1D_t  q2Lim_W (void) const;  ///< q2 limits @ fixed W
  Range1D_t  Q2Lim   (void) const;  ///< Q2 limits
  Range1D_t  q2Lim   (void) const;  ///< q2 limits
  Range1D_t  XLim    (void) const;  ///< x  limits
  Range1D_t  YLim    (void) const;  ///< y  limits
  Range1D_t  YLim_X  (void) const;  ///< y  limits @ fixed x
  Range1D_t  YLim    (double xsi) const;  ///< y  limits (COH)
  Range1D_t  YLim_X  (double xsi) const;  ///< y  limits @ fixed x (COH)
  Range1D_t  TLim    (void) const;  ///< t  limits
  Range1D_t  WLim_SPP(void) const;          ///< W  limits for single pion production models
  Range1D_t  WLim_SPP_iso (void) const;     ///< W  limits for resonance single pion production on isoscalar nucleon
  Range1D_t  Q2Lim_W_SPP  (void) const;     ///< Q2 limits @ fixed W for single pion production models
  Range1D_t  Q2Lim_W_SPP_iso (void) const;  ///< Q2 limits @ fixed W for resonance single pion production on isoscalar nucleon

  static double GetTMaxDFR();

private:
  void Init(void);

  const Interaction * fInteraction;

ClassDef(KPhaseSpace,2)
};

}      // genie namespace
#endif // _KINE_PHASE_SPACE_H_
