//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 26, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Interaction/IUtils.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
Range1D_t  genie::utils::kinematics::KineRange(
                                   const Interaction * const i, KineVar_t k)
{
  switch(k) {
  case(kKVW)  : return WRange(i);  break;
  case(kKVQ2) : return Q2Range(i); break;
  case(kKVq2) : return q2Range(i); break;
  default:
   SLOG("KineLimits", pERROR) 
              << "Couldn't compute limits for " 
                                << KineVar::AsString(k) << "using \n"<< *i;
   Range1D_t R(-1.,-1);
   return R;
  }
}
//____________________________________________________________________________
double genie::utils::kinematics::PhaseSpaceVolume(
                          const Interaction * const in, KinePhaseSpace_t ps)
{
  double vol = 0;

  switch(ps) {

  case(kPSQ2fE):
  {
    Range1D_t Q2 = KineRange(in, kKVQ2);
    vol = Q2.max - Q2.min;
    return vol;
    break;
  }
  case(kPSq2fE):
  {
    Range1D_t q2 = KineRange(in, kKVq2);
    vol = q2.min - q2.max;
    return vol;
    break;
  }
  case(kPSWfE):
  {
    Range1D_t W = KineRange(in, kKVW);
    vol = W.max - W.min;
    return vol;
    break;
  }
  case(kPSWQ2fE):
  {
    Range1D_t W = KineRange(in, kKVW);
    if(W.max<0) return 0;
    const int kNW = 100;  
    double dW = (W.max-W.min)/(kNW-1);	
    Interaction interaction(*in);
    for(int iw=0; iw<kNW; iw++) {
      interaction.GetKinematicsPtr()->SetW(W.min + iw*dW);
      Range1D_t Q2 = KineRange(&interaction, kKVQ2);
      double dQ2 = (Q2.max-Q2.min);	
      vol += (dW*dQ2);
    }
    return vol;
    break;
  }
  case(kPSxyfE):
  {
    Range1D_t W = KineRange(in, kKVW);
    if(W.max<0) return 0;

    const InitialState & init_state = in->GetInitialState();
    double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
    double M  = init_state.GetTarget().StruckNucleonP4()->M(); 

    const int    kNx = 100;  
    const int    kNy = 100;  
    const double kdx = (kMaxX-kMinX)/(kNx-1);	
    const double kdy = (kMaxY-kMinY)/(kNy-1);	
    const double kdV = kdx*kdy;

    double cW=-1, cQ2 = -1;
   
    Interaction interaction(*in);

    for(int ix=0; ix<kNx; ix++) {
      double x = kMinX+ix*kdx;
      for(int iy=0; iy<kNy; iy++) {
         double y = kMinY+iy*kdy;

         XYtoWQ2(Ev, M, cW, cQ2, x, y);
         if(!math::IsWithinLimits(cW, W)) continue;

         interaction.GetKinematicsPtr()->SetW(cW);
         Range1D_t Q2 = KineRange(&interaction, kKVQ2);
         if(!math::IsWithinLimits(cQ2, Q2)) continue;

         vol += kdV;
      }
    }
    return vol;
    break;
  }
  default:
   SLOG("KineLimits", pERROR) 
     << "Couldn't compute phase space volume for " 
                    << KinePhaseSpace::AsString(ps) << "using \n" << *in;
   return 0;
  }
}
//____________________________________________________________________________
double genie::utils::kinematics::Jacobian(
            const Interaction * const i, 
                              KinePhaseSpace_t fromps, KinePhaseSpace_t tops)
{
// returns the Jacobian for the transformation: fromps -> tops
//
  SLOG("KineLimits", pDEBUG) 
       << "Computing Jacobian for transformation: "
                << KinePhaseSpace::AsString(fromps) << " --> " 
                                            << KinePhaseSpace::AsString(tops);
  double J=0;
  bool forward;
  const Kinematics & kine = i->GetKinematics();

  // ****** transformation: {Q2}|E -> {lnQ2}|E
  if ( TransformMatched(fromps,tops,kPSQ2fE,kPSlogQ2fE,forward) )
  {
    J = kine.Q2();
  } 
  // ****** transformation: {Q2}|E -> {QD2}|E
  else if ( TransformMatched(fromps,tops,kPSQD2fE,kPSQ2fE,forward) ) 
  {
    J = TMath::Power(1+kine.Q2()/kMQD2,-2)/kMQD2;
  } 
  // ****** transformation: {W,Q2}|E -> {W,QD2}|E
  else if ( TransformMatched(fromps,tops,kPSWQD2fE,kPSWQ2fE,forward) ) 
  {
    J = TMath::Power(1+kine.Q2()/kMQD2,-2)/kMQD2;
  } 
  // ****** transformation: {x,y}|E -> {lnx,lny}|E
  else if ( TransformMatched(fromps,tops,kPSxyfE,kPSlogxlogyfE,forward) ) 
  {
    J = kine.x() * kine.y();
  } 
  // ****** transformation: {W,Q2}|E -> {W,lnQ2}|E
  else if ( TransformMatched(fromps,tops,kPSWQ2fE,kPSWlogQ2fE,forward) ) 
  {
    J = kine.Q2();
  } 
  // ****** transformation: {W2,Q2}|E -> {x,y}|E
  else if ( TransformMatched(fromps,tops,kPSW2Q2fE,kPSxyfE,forward) ) 
  {
    const InitialState & init_state = i->GetInitialState();
    double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
    double M  = init_state.GetTarget().StruckNucleonP4()->M(); 
    double y  = kine.y();
    J = TMath::Power(2*M*Ev,2) * y;
  } 
  // ****** transformation: {W2,logQ2}|E -> {x,y}|E
  else if ( TransformMatched(fromps,tops,kPSW2logQ2fE,kPSxyfE,forward) ) 
  {
    J = Jacobian(i,kPSW2Q2fE,kPSxyfE)/kine.Q2();
  } 
  else {
     SLOG("KineLimits", pFATAL) 
       << "*** Can not compute Jacobian for transforming: "
                     << KinePhaseSpace::AsString(fromps) << " <--> " 
                                            << KinePhaseSpace::AsString(tops);
     exit(1);
  }

  // if any of the above transforms was reverse, invert the jacobian
  if(!forward) J = 1./J;

  return J;
}
//____________________________________________________________________________
bool genie::utils::kinematics::TransformMatched(
	       KinePhaseSpace_t inpa, KinePhaseSpace_t inpb,
                           KinePhaseSpace_t a, KinePhaseSpace_t b, bool & fwd)
{
  
  bool matched = (a==inpa&&b==inpb || a==inpb&&b==inpa); // match a->b or b->a
  if(matched) {
    if(a==inpa&&b==inpb) fwd = true;  // matched forward transform: a->b
    else                 fwd = false; // matched reverse transform: b->a
  }
  return matched;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::WRange(
                                        const Interaction * const interaction)
{
// Computes kinematical limits for W in DIS and RES v interactions
//
  Range1D_t W;

  const InitialState & init_state = interaction->GetInitialState();

  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonP4()->M(); //can be off m/shell
  double ml = interaction->GetFSPrimaryLepton()->Mass();
  double M2 = M*M;
  double s  = M2 + 2*M*Ev;

  assert (s>0);

  W.min  = kNeutronMass + kPionMass;
  W.max  = TMath::Sqrt(s)-ml;

  if(W.max<=W.min) {
    W.min = -1;
    W.max = -1;
  }

  return W;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range(
                                       const Interaction * const interaction)
{
// Computes kinematical limits for Q2 (>0)

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  const ProcessInfo & pi = interaction->GetProcessInfo();

  bool handle = pi.IsDeepInelastic() || pi.IsResonant() || pi.IsQuasiElastic();

  if(!handle) return Q2;

  const InitialState & init_state  = interaction->GetInitialState();

  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);
  double M   = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m/shell
  double ml  = interaction->GetFSPrimaryLepton()->Mass();
  double M2  = M*M;
  double ml2 = ml*ml;
  double s   = M2 + 2*M*Ev;

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;

  assert (s>0);

  double W  = CalcW(interaction);
  double W2 = W*W;

  if(pi.IsDeepInelastic() || pi.IsResonant()) {

    Range1D_t rW = WRange(interaction);
    SLOG("KineLimits", pDEBUG)
                 << "W  = " << W << ", [Wmin = " << rW.min
                                            << ", Wmax = " << rW.max << "]";
    if(W < rW.min || W > rW.max) return Q2;
  }

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
  if(Q2.max < Q2.min) {Q2.min = -1; Q2.max = -1;}

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::q2Range(
                                      const Interaction * const interaction)
{
// Computes kinematical limits for q2 (<0) 

  Range1D_t Q2 = utils::kinematics::Q2Range(interaction);

  Range1D_t q2;

  q2.min = - Q2.max;
  q2.max = - Q2.min;

  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Q2Range_W(
                          const Interaction * const interaction, Range1D_t rW)
{
  Range1D_t Q2;
  Q2.min = 999;
  Q2.max = -999;

  unsigned int N=100;
  double dW = (rW.max-rW.min) / (double)(N-1);

  for(unsigned int i=0; i<N; i++) {
    double W = rW.min + i*dW;
    interaction->GetKinematicsPtr()->SetW(W);
    Range1D_t rQ2 = Q2Range(interaction);

    if(rQ2.min>=0) Q2.min = TMath::Min(Q2.min, rQ2.min);
    if(rQ2.max>=0) Q2.max = TMath::Max(Q2.max, rQ2.max);
  }
  return Q2;
}
//____________________________________________________________________________
void genie::utils::kinematics::MinXY(
                        const Interaction * const in, double & x, double & y)
{
// Compute the minimum x,y for the current interaction (calculated as the x,y
// at the minimum W,Q2)

  Interaction interaction(*in);

  Range1D_t W = KineRange(&interaction,kKVW);
  interaction.GetKinematicsPtr()->SetW(W.min);

  Range1D_t Q2 = KineRange(&interaction,kKVQ2);

  const InitialState & init_state =  interaction.GetInitialState();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonP4()->M(); //can be off m/shell

  WQ2toXY(Ev,M,W.min,Q2.min,x,y);
}
//____________________________________________________________________________
double genie::utils::kinematics::EnergyThreshold(
                                        const Interaction * const interaction)
{
  const Target &      tgt  = interaction->GetInitialState().GetTarget();
  const ProcessInfo & proc = interaction->GetProcessInfo();

  double ml = interaction->GetFSPrimaryLepton()->Mass();

  if (proc.IsCoherent()) {

    int tgtpdgc = tgt.PDGCode(); // nuclear target PDG code (1AAAZZZ000)
    double MA   = PDGLibrary::Instance()->Find(tgtpdgc)->Mass(); 
    double m    = ml + kPionMass;
    double m2   = TMath::Power(m,2);
    double Ethr = m + 0.5*m2/MA;
    return Ethr;
  }

  if(proc.IsQuasiElastic() || proc.IsResonant() || proc.IsDeepInelastic()) {

    assert( tgt.StruckNucleonIsSet() );
    double Mn   = tgt.StruckNucleonP4()->M();
    double Mn2  = TMath::Power(Mn,2);
    double Wmin = (proc.IsQuasiElastic()) ? Mn : kNeutronMass+kPionMass;
    double smin = TMath::Power(Wmin+ml,2.);
    double Ethr = 0.5*(smin-Mn2)/Mn;
    return Ethr;
  }

  SLOG("KineLimits", pERROR)
      	      << "Can't compute threshold for \n" << *interaction;
  exit(1);

  return 99999999;
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

  double Mn   = init_state.GetTarget().StruckNucleonP4()->M();
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
double genie::utils::kinematics::Q2toQD2(double Q2)
{
// Q2 -> QD2 transformation. QD2 takes out the dipole form of the form factors
// making the differential cross section to be flatter and speeding up the
// kinematical selection.

  assert(Q2>0);
  return TMath::Power(1+Q2/kMQD2, -1);
}
//____________________________________________________________________________
double genie::utils::kinematics::QD2toQ2(double QD2)
{
  assert(QD2>0);
  return kMQD2*(1/QD2-1);
}
//____________________________________________________________________________
double genie::utils::kinematics::CalcQ2(const Interaction * const interaction)
{
  const Kinematics & kinematics = interaction->GetKinematics();

  if (kinematics.KVSet(kKVQ2) || kinematics.KVSet(kKVq2)) {
    double Q2 = kinematics.Q2();
    return Q2;
  }
  if (kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->GetInitialState();
    double Mn = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m/shell
    double x  = kinematics.x();
    double y  = kinematics.y();
    double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
    double Q2 = 2*Mn*Ev*x*y;
    return Q2;
  }

  SLOG("KineLimits", pERROR) << "Couldn't compute Q^2 for \n"<< *interaction;
  return 0;
}
//____________________________________________________________________________
double genie::utils::kinematics::CalcW(const Interaction * const interaction)
{
  const ProcessInfo & process_info = interaction->GetProcessInfo();

  if(process_info.IsQuasiElastic()) {
    //    const InitialState & init_state  = interaction->GetInitialState();
    //    double M  = init_state.GetTarget().StruckNucleonP4()->M(); 
    //    return M;

    // hadronic inv. mass is equal to the recoil nucleon on-shell mass
    int    rpdgc = utils::interaction::RecoilNucleonPdgCode(interaction);
    double M = PDGLibrary::Instance()->Find(rpdgc)->Mass();
    return M;
  }

  const Kinematics & kinematics = interaction->GetKinematics();

  if(kinematics.KVSet(kKVW)) {
    double W = kinematics.W();
    return W;
  }

  if(kinematics.KVSet(kKVx) && kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->GetInitialState();
    double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
    double M  = init_state.GetTarget().StruckNucleonP4()->M(); 
    double M2 = M*M;
    double x  = kinematics.x();
    double y  = kinematics.y();
    double W2 = TMath::Max(0., M2 + 2*Ev*M*y*(1-x));
    double W  = TMath::Sqrt(W2);
    return W;
  }

  SLOG("KineLimits", pERROR) << "Couldn't compute W for \n"<< *interaction;
  return 0;
}
//____________________________________________________________________________
void genie::utils::kinematics::WQ2toXY(
            double Ev, double M, double W, double Q2, double & x, double & y)
{
// Converts (W,Q2) => (x,y)
// Uses the system: a) W^2 - M^2 = 2*Ev*M*y*(1-x) and b) Q^2 = 2*x*y*M*Ev
// Ev is the neutrino energy at the struck nucleon rest frame
// M is the nucleon mass - it does not need to be on the mass shell

  double M2  = TMath::Power(M,2);
  double W2  = TMath::Power(W,2);

  x = Q2 / (W2-M2+Q2);
  y = (W2-M2+Q2) / (2*M*Ev);

  x = TMath::Min(1.,x);
  y = TMath::Min(1.,y);
  x = TMath::Max(0.,x);
  y = TMath::Max(0.,y);

  LOG("KineLimits", pDEBUG) 
        << "(W=" << W << ",Q2=" << Q2 << ") => (x="<< x << ", y=" << y<< ")";
}
//___________________________________________________________________________
void genie::utils::kinematics::XYtoWQ2(
            double Ev, double M, double & W, double & Q2, double x, double y)
{
// Converts (x,y) => (W,Q2)
// Uses the system: a) W^2 - M^2 = 2*Ev*M*y*(1-x) and b) Q^2 = 2*x*y*M*Ev
// Ev is the neutrino energy at the struck nucleon rest frame
// M is the nucleon mass - it does not need to be on the mass shell

  double M2  = TMath::Power(M,2);
  double W2  = M2 + 2*Ev*M*y*(1-x);

  W  = TMath::Sqrt(TMath::Max(0., W2));
  Q2 = 2*x*y*M*Ev;

  LOG("KineLimits", pDEBUG) 
      << "(x=" << x << ",y=" << y << " => (W=" << W << ",Q2=" << Q2 << ")";
}
//___________________________________________________________________________
double genie::utils::kinematics::XYtoW(
                                     double Ev, double M, double x, double y)
{
// Converts (x,y) => W
// Ev is the neutrino energy at the struck nucleon rest frame
// M is the nucleon mass - it does not need to be on the mass shell

  double M2  = TMath::Power(M,2);
  double W2  = M2 + 2*Ev*M*y*(1-x);
  double W  = TMath::Sqrt(TMath::Max(0., W2));

  LOG("KineLimits", pDEBUG) << "(x=" << x << ",y=" << y << ") => W=" << W;

  return W;
}
//___________________________________________________________________________
double genie::utils::kinematics::XYtoQ2(
                                     double Ev, double M, double x, double y)
{
// Converts (x,y) => Q2
// Ev is the neutrino energy at the struck nucleon rest frame
// M is the nucleon mass - it does not need to be on the mass shell

  double Q2 = 2*x*y*M*Ev;

  LOG("KineLimits", pDEBUG) << "(x=" << x << ",y=" << y << ") => Q2=" << Q2;

  return Q2;
}
//___________________________________________________________________________
double genie::utils::kinematics::SlowRescalingVar(
                             const Interaction * const interaction, double mc)
{
  const InitialState & init_state = interaction -> GetInitialState();
  const Kinematics &   kinematics = interaction -> GetKinematics();

  double mc2 = TMath::Power(mc,2);
  double Mn  = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m/shell
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
double genie::utils::kinematics::RESImportanceSamplingEnvelope(
                                                     double * x, double * par)
{
  //-- inputs
  double QD2   = x[0];   // QD2 (Q2 transformed to take out the dipole form)
  double W     = x[1];   // invariant mass

  //-- parameters
  double mD    = par[0]; // resonance mass
  double gD    = par[1]; // resonance width
  double xsmax = par[2]; // safety factor * max cross section in (W,Q2)
  double Wmax  = par[3]; // kinematically allowed Wmax

  double func = 0;

  if(Wmax > mD) {
     // -- if the resonance mass is within the kinematical range then
     // -- the envelope is defined as a plateau above the resonance peak,
     // -- a steeply falling leading edge (low-W side) and a more slowly
     // -- falling trailing edge (high-W side)
     if(W < mD-gD/2) {
       //low-W falling edge
       double plateau_edge = mD-gD/2;
       func = xsmax / (1 + 5* TMath::Power((W-plateau_edge)/gD,2));
     } else if (W > mD+gD/2) {
       //high-W falling edge
       double plateau_edge = mD+gD/2;
       func = xsmax / (1 + 2 * TMath::Power((W-plateau_edge)/gD,2));
     } else {
       // plateau
       func = xsmax;
     }
  } else {
     // -- if the resonance mass is above the kinematical range then the
     // -- envelope is a small plateau just bellow Wmax and a falling edge
     // -- at lower W
     if (W > Wmax-0.1) {
       // plateau
       func = xsmax;
     } else {
       //low-W falling edge
       double plateau_edge = Wmax-0.1;
       func = xsmax / (1 + TMath::Power((W-plateau_edge)/gD,2));
     } 
  }

  if(QD2<0.5) {
    func *= (1 - (0.5-QD2));
  }

  return func;
}
//___________________________________________________________________________
double genie::utils::kinematics::DISImportanceSamplingEnvelope(
                                                     double * x, double * par)
{
  //-- inputs
  double xb = x[0];   // scaling variable x brorken
  //double y  = x[1];   // inelasticity y

  //-- parameters
  double xpeak = par[0]; // x peak
  double xmin  = par[1]; // min x
  double xmax  = 1.;     // max x
  double xsmax = par[2]; // safety factor * max cross section in (x,y)

  double func = 0;

  if(xb < xpeak/20.) {
       //low-x falling edge
       double plateau_edge = xpeak/20.;
       double slope = xsmax/(xmin - plateau_edge);
       func = xsmax - slope * (xb - plateau_edge);
  } else if (xb > 2*xpeak) {
       //high-x falling edge
       double plateau_edge = 2*xpeak;
       double slope = xsmax/(xmax - plateau_edge);
       func = xsmax - slope * (xb - plateau_edge);
  } else {
       // plateau
       func = xsmax;
  } 
  return func;
}
//___________________________________________________________________________
double genie::utils::kinematics::COHImportanceSamplingEnvelope(
                                                    double * x, double * par)
{
  //-- inputs
  double xb = x[0];   // x
  double y  = x[1];   // y

  //-- parameters
  double ymin  = par[0]; // min y
  double ymax  = par[1]; // max y
  double xsmax = par[2]; // safety factor * max cross section in (x,y)
  double Ev    = par[3]; // neutrino energy;

  if(y<=ymin || y>=ymax) return 0.;
  if(xb<=0.  || xb>=1.)  return 0.;

  double EpiP = 0.26;
  double yP   = EpiP/Ev;

  double func = 0;
  if(y>yP) {
    // falling edge
    double d = TMath::Power( (y-yP)/0.015, 2);
    func = xsmax / (1 + d);
  } else {
    // plateau
    func = xsmax;
  } 
  return func;
}
//___________________________________________________________________________


