//____________________________________________________________________________
/*!

\namespace  genie::utils::kinematics

\brief      Kinematical utilities

\author     Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
            University of Liverpool & STFC Rutherford Appleton Laboratory

            Changes required to implement the GENIE Boosted Dark Matter module
            were installed by Josh Berger (Univ. of Wisconsin)

\created    November 26, 2004

\cpright    Copyright (c) 2003-2023, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org            
*/
//____________________________________________________________________________

#ifndef _KINE_UTILS_H_
#define _KINE_UTILS_H_

#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Utils/Range1.h"

namespace genie {
namespace utils {

namespace kinematics
{
  //-- methods used for computing phase space volumes
  double PhaseSpaceVolume (const Interaction * const i, KinePhaseSpace_t ps);

  //-- methods used for computing the Jacobians associated with phase space transformations
  double Jacobian (const Interaction * const i, KinePhaseSpace_t f, KinePhaseSpace_t t);
  bool   TransformMatched (KinePhaseSpace_t ia, KinePhaseSpace_t ib,
                                    KinePhaseSpace_t a, KinePhaseSpace_t b, bool & fwd);

  //-- kinematical limits
  Range1D_t  InelWLim    (double Ev, double M, double ml);
  Range1D_t  InelQ2Lim_W (double Ev, double M, double ml, double W, double Q2min_cut =    controls::kMinQ2Limit);
  Range1D_t  Inelq2Lim_W (double Ev, double M, double ml, double W, double q2min_cut = -1*controls::kMinQ2Limit);
  Range1D_t  InelQ2Lim   (double Ev, double M, double ml, double Q2min_cut =    controls::kMinQ2Limit);
  Range1D_t  Inelq2Lim   (double Ev, double M, double ml, double q2min_cut = -1*controls::kMinQ2Limit);
  Range1D_t  InelXLim    (double Ev, double M, double ml);
  Range1D_t  InelYLim    (double Ev, double M, double ml);
  Range1D_t  InelYLim_X  (double Ev, double M, double ml, double x);
  Range1D_t  CohW2Lim    (double Mn, double m_produced, double mlep, double Ev, double Q2);
  Range1D_t  CohNuLim    (double W2min, double W2max, double Q2, double Mn, double xsi);
  Range1D_t  CohYLim     (double Mn, double m_produced, double mlep, double Ev, double Q2, double xsi);
  Range1D_t  CohYLim     (double EvL, double ml);
  Range1D_t  CohXLim     (void);
  Range1D_t  CohQ2Lim    (double Mn, double m_produced, double mlep, double Ev);
  Range1D_t  Cohq2Lim    (double Mn, double m_produced, double mlep, double Ev);
  Range1D_t  CEvNSQ2Lim  (double Ev);
  Range1D_t  DarkWLim    (double Ev, double M, double ml);
  Range1D_t  DarkQ2Lim_W (double Ev, double M, double ml, double W, double Q2min_cut =    controls::kMinQ2Limit);
  Range1D_t  Darkq2Lim_W (double Ev, double M, double ml, double W, double q2min_cut = -1*controls::kMinQ2Limit);
  Range1D_t  DarkQ2Lim   (double Ev, double M, double ml, double Q2min_cut =    controls::kMinQ2Limit);
  Range1D_t  Darkq2Lim   (double Ev, double M, double ml, double q2min_cut = -1*controls::kMinQ2Limit);
  Range1D_t  DarkXLim    (double Ev, double M, double ml);
  Range1D_t  DarkYLim    (double Ev, double M, double ml);
  Range1D_t  DarkYLim_X  (double Ev, double M, double ml, double x);

  //-- helpers for kinematic limits
  double CohW2Min(double Mn, double m_produced);

  //-- kinematical variable transforms
  double QD2toQ2 (double QD2);
  double Q2toQD2 (double Q2);
  void   WQ2toXY (double Ev, double M, double W, double Q2, double & x, double & y);
  void   XYtoWQ2 (double Ev, double M, double & W, double & Q2, double x, double y);
  void   XQ2toWY (double Ev, double M, double & W, double Q2, double x, double & y);
  double XYtoW   (double Ev, double M, double x, double y);
  double XYtoQ2  (double Ev, double M, double x, double y);
  double Q2YtoX  (double Ev, double M, double Q2, double y);

  void  UpdateWQ2FromXY(const Interaction * in);
  void  UpdateXYFromWQ2(const Interaction * in);
  void  UpdateWYFromXQ2(const Interaction * in);
  void  UpdateXFromQ2Y(const Interaction * in);

  //-- methods used to apply cuts to kinematical limits
  void ApplyCutsToKineLimits (Range1D_t & r, double min, double max);

  double Q2 (const Interaction * const i);
  double W  (const Interaction * const i);

  //-- charm threshold and slow rescaling variable
  bool   IsAboveCharmThreshold (double x, double Q2, double M, double mc);
  double SlowRescalingVar      (double x, double Q2, double M, double mc);

  //-- Functions used to define differential cross section distribution envelopes 
  //-- for importance sampling in kinematical selection modules
  double RESImportanceSamplingEnvelope(double * x, double * par);
  double DISImportanceSamplingEnvelope(double * x, double * par);
  double COHImportanceSamplingEnvelope(double * x, double * par);

  namespace electromagnetic
  {
   Range1D_t  InelWLim    (double El, double ml, double M);
   Range1D_t  InelQ2Lim_W (double El, double ml, double M, double W);
   Range1D_t  Inelq2Lim_W (double El, double ml, double M, double W);
   Range1D_t  InelQ2Lim   (double El, double ml, double M);
   Range1D_t  Inelq2Lim   (double El, double ml, double M);
   Range1D_t  InelXLim    (double El, double ml, double M);
   Range1D_t  InelYLim    (double El, double ml, double M);
   Range1D_t  InelYLim_X  (double El, double ml, double M, double x);

   static const double kMinQ2Limit   = 0.02;  // GeV^2 // Q2 threshold relevant for em scattering events
  }

} // kinematics namespace
} // utils namespace
} // genie namespace

#endif // _KINE_UTILS_H_
