//____________________________________________________________________________
/*!

\namespace  genie::utils::kinematics

\brief      Kinematical limits for DIS,QEL,RES

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    Novemner 26, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Interaction/KineVar.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
Range1D_t genie::utils::kinematics::WRange(const Interaction * const interaction)
{
// Computes kinematical limits for W in DIS and RES v interactions
//
  Range1D_t W;

  const InitialState & init_state = interaction->GetInitialState();

  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonMass();
  double ml = interaction->GetFSPrimaryLepton()->Mass();
  double M2 = M*M;
  double s  = M2 + 2*M*Ev;

  assert (s>0);

  W.min  = kNeutronMass + kPionMass;
  W.max  = TMath::Sqrt(s)-ml;

  return W;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range(const Interaction * const interaction)
{
// Computes kinematical limits for Q2 (>0) in QEL, DIS or RES v interactions.

  Range1D_t Q2;

  const ProcessInfo & process_info = interaction->GetProcessInfo();

  if ( process_info.IsQuasiElastic()  ) {

       Q2 = utils::kinematics::Q2Range_M(interaction);

  } else if ( process_info.IsDeepInelastic() ) {

       Q2 = utils::kinematics::Q2Range_xy(interaction);

  } else if ( process_info.IsResonant() ) {

       Q2 = utils::kinematics::Q2Range_W(interaction);

  } else {

       Q2.min = 0;
       Q2.max = 0;
  }

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::q2Range(const Interaction * const interaction)
{
// Computes kinematical limits for q2 (<0) in QEL, DIS or RES v interactions.

  Range1D_t Q2 = utils::kinematics::Q2Range(interaction);

  Range1D_t q2;

  q2.min = - Q2.max;
  q2.max = - Q2.min;

  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range_W(const Interaction * const interaction)
{
// Computes kinematical limits for Q2 (>0) in [RES] v interactions.
// The Q^2 range depends on W which is extracted directly from the input
// interaction

  SLOG("KineLimits", pDEBUG) << *interaction;

  Range1D_t Q2;

  Q2.min = -1;
  Q2.max = -1;

  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev    = p4->Energy();
  double M     = init_state.GetTarget().StruckNucleonMass();
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double M2    = M*M;
  double ml2   = ml*ml;
  double s     = M2 + 2*M*Ev;
  double W     = interaction->GetKinematics().W();
  double W2    = W*W;

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;

  assert (s>0);

  Range1D_t rW = WRange(interaction);

  SLOG("KineLimits", pDEBUG)
                << "W  = " << W << ", [Wmin = " << rW.min
                                         << ", Wmax = " << rW.max << "]";
  if( W >= rW.min && W <= rW.max) {

     double auxC = 0.5*(s-M2)/s;

     double aux1 = s + ml2 - W2;
     double aux2 = aux1*aux1 - 4*s*ml2;

     (aux2 < 0) ? ( aux2 = 0 ) : ( aux2 = TMath::Sqrt(aux2) );

     Q2.max = -ml2 + auxC * (aux1 + aux2); // => 0
     Q2.min = -ml2 + auxC * (aux1 - aux2); // => 0

     // guard against overflows
     Q2.max = TMath::Max(0., Q2.max);
     Q2.min = TMath::Max(0., Q2.min);

     // limit the minimum Q2
     if(Q2.min < kMinQ2Limit) Q2.min = kMinQ2Limit;
     if(Q2.max < kMinQ2Limit) Q2.max = kMinQ2Limit;
  }

  delete p4;

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range_xy(
                                        const Interaction * const interaction)
{
// Computes kinematical limits for Q2 (>0) in [DIS] v interactions
// They depend on the W which is computed from the x and y kinematical vars
// that are extracted from the input interaction
//
  SLOG("KineLimits", pDEBUG) << *interaction;

  Range1D_t Q2;

  Q2.min = -1;
  Q2.max = -1;

  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev    = p4->Energy();

  double x     = interaction->GetKinematics().x();
  double y     = interaction->GetKinematics().y();

  double M     = init_state.GetTarget().StruckNucleonMass();
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double M2    = M*M;
  double ml2   = ml*ml;
  double s     = M2 + 2*M*Ev;

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;

  assert (s>0);

  double Wmin  = kNeutronMass + kPionMass;
  double Wmax  = TMath::Sqrt(s)-ml;

  double W2 = TMath::Max(0., M2 + 2*Ev*M*y*(1-x));
  double W  = TMath::Sqrt(W2);

  SLOG("KineLimits", pDEBUG)
                << "W  = " << W << ", [Wmin = " << Wmin
                                               << ", Wmax = " << Wmax << "]";
  if(W >= Wmin & W <= Wmax) {

     double tmp1  = s+ml2-W2;
     double tmp2  = TMath::Sqrt( TMath::Max(0., tmp1*tmp1 - 4*s*ml2) );

     Q2.min = - ml2 + 0.5 * ((s-M2)/s) * (tmp1-tmp2);
     Q2.max = - ml2 + 0.5 * ((s-M2)/s) * (tmp1+tmp2);

     // guard against overflows
     Q2.max = TMath::Max(0., Q2.max);
     Q2.min = TMath::Max(0., Q2.min);

     // limit the minimum Q2
     if(Q2.min < kMinQ2Limit) Q2.min = kMinQ2Limit;
     if(Q2.max < kMinQ2Limit) Q2.max = kMinQ2Limit;
  }

  delete p4;

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range_M(const Interaction * const interaction)
{
// Computes kinematical limits for Q2 (>0) in [QEL] v interactions
// For QEL events W = M
//
  SLOG("KineLimits", pDEBUG) << *interaction;

  Range1D_t Q2;

  Q2.min = -1;
  Q2.max = -1;

  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Ev    = p4->Energy();

  double M     = init_state.GetTarget().StruckNucleonMass();
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double M2    = M*M;
  double ml2   = ml*ml;
  double s     = M2 + 2*M*Ev;

  assert (s>0);

  double tmp1  = s+ml2-M2;
  double tmp2  = TMath::Sqrt( TMath::Max(0., tmp1*tmp1 - 4*s*ml2) );

  Q2.min = - ml2 + 0.5 * ((s-M2)/s) * (tmp1-tmp2);
  Q2.max = - ml2 + 0.5 * ((s-M2)/s) * (tmp1+tmp2);

  // guard against overflows
  Q2.max = TMath::Max(0., Q2.max);
  Q2.min = TMath::Max(0., Q2.min);

  // limit the minimum Q2
  if(Q2.min < kMinQ2Limit) Q2.min = kMinQ2Limit;
  if(Q2.max < kMinQ2Limit) Q2.max = kMinQ2Limit;

  delete p4;

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::q2Range_W(const Interaction * const interaction)
{
// Computes kinematical limits for q2 (<0) in [RES] v interactions.
// The q^2 range depends on W which is extracted directly from the input
// interaction

  Range1D_t Q2 = utils::kinematics::Q2Range_W(interaction);

  Range1D_t q2;

  q2.min = - Q2.max;
  q2.max = - Q2.min;

  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::q2Range_xy(
                                        const Interaction * const interaction)
{
// Computes kinematical limits for q2 (<0) in [DIS] v interactions
// They depend on the W which is computed from the x and y kinematic variables
// that are extracted from the input interaction
//

  Range1D_t Q2 = utils::kinematics::Q2Range_xy(interaction);

  Range1D_t q2;

  q2.min = - Q2.max;
  q2.max = - Q2.min;

  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::q2Range_M(const Interaction * const interaction)
{
// Computes kinematical limits for q2 (<0) in [QEL] v interactions
// For QEL events W = M
//
  Range1D_t Q2 = utils::kinematics::Q2Range_M(interaction);

  Range1D_t q2;

  q2.min = - Q2.max;
  q2.max = - Q2.min;

  return q2;
}
//____________________________________________________________________________
double genie::utils::kinematics::EnergyThreshold(
                                        const Interaction * const interaction)
{
  double Ethr  = 0;
  double ml    = interaction->GetFSPrimaryLepton()->Mass();
  double ml2   = ml*ml;
  double Mnuc  = kNucleonMass;
  double Mnuc2 = Mnuc*Mnuc;

  const ProcessInfo & proc = interaction->GetProcessInfo();

  if(proc.IsResonant() || proc.IsDeepInelactic()) {
      Range1D_t rW    = utils::kinematics::WRange(interaction);
      double    Wmin  = rW.min;
      double    Wmin2 = Wmin*Wmin;

      Ethr = (ml2 + Wmin2 + 2*ml*Wmin - Mnuc2) / (2*Mnuc);
  }
  else if (proc.IsQuasiElastic()) {
      Ethr = (ml2 + 2*ml*Mnuc) / (2*Mnuc);
  }
  else {
      SLOG("KineLimits", pERROR)
             << "Doesn't compute Ethreshold for this interaction":
  }

  return Ethreshold;
}
//____________________________________________________________________________
bool genie::utils::kinematics::IsAboveCharmThreshold(
                             const Interaction * const interaction, double mc)
{
// check along the lines of the NeuGEN check in nu_dissf
// mc is the input charm mass - make sure that a mc scale is used consistently
// throughtout GENIE
//
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  int lightest_charm_hadron = kPdgDMinus; // c=+/-1, s=0

  double Mn   = init_state.GetTarget().StruckNucleonMass();
  double Mn2  = TMath::Power(Mn,2);
  double Q2   = utils::kinematics::CalcQ2(interaction);
  double x    = kinematics.x();
  double v    = 0.5*Q2/(Mn*x);
  double W2   = TMath::Max(0., Mn2+2*Mn*v-Q2);
  double W    = TMath::Sqrt(W2);
  double mD   = PDGLibrary::Instance()->Find(lightest_charm_hadron)->Mass();
  double Wmin = Mn+mD;
  double xc   = utils::kinematics::SlowRescalingVar(interaction,mc);

  if(xc>=1 || W<=Wmin) return false;
  else                 return true;
}
//____________________________________________________________________________
double genie::utils::kinematics::CalcQ2(const Interaction * const interaction)
{
  double Q2=0.;
  const Kinematics & kinematics = interaction->GetKinematics();
  if( kinematics.KVSet(kKVQ2) || kinematics.KVSet(kKVq2)) Q2 = kinematics.Q2();
  else if ( kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->GetInitialState();
    double Mn = init_state.GetTarget().StruckNucleonMass();
    double x  = kinematics.x();
    double y  = kinematics.y();
    double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
    Q2 = 2*Mn*Ev*x*y;
  }
  return Q2;
}
//____________________________________________________________________________
double genie::utils::kinematics::SlowRescalingVar(
                             const Interaction * const interaction, double mc)
{
  const InitialState & init_state = interaction -> GetInitialState();
  const Kinematics &   kinematics = interaction -> GetKinematics();

  double mc2 = TMath::Power(mc,2);
  double Mn  = init_state.GetTarget().StruckNucleonMass();
  double x   = kinematics.x();
  double Q2  = utils::kinematics::CalcQ2(interaction);
  double v   = 0.5*Q2/(Mn*x);
  double xc  = x + 0.5*mc2/(Mn*v);

  return xc;
}
//____________________________________________________________________________
void genie::utils::kinematics::ApplyCutsToKineLimits(
                         Range1D_t & range, double min_cut, double max_cut)
{
  // if the min,max are within the existing limits, the cut can be applied
  // by narrowing down the xisting limits
  if ( utils::math::IsWithinLimits(min_cut, range ) ) range.min = min_cut;
  if ( utils::math::IsWithinLimits(max_cut, range ) ) range.max = max_cut;

  // if the min-cut is above the existing max-limit or
  // if the max-cut is below the existing min-limit then
  // the range should be invalidated

  if (min_cut > range.max || max_cut < range.min) {

    range.min = 0;
    range.max = 0;
  }
}
//____________________________________________________________________________

