//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

	  Afroditi Papadopoulou <apapadop \at mit.edu>
	  Massachusetts Institute of Technology

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <cstdlib>
#include <limits>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::kinematics::PhaseSpaceVolume(
  const Interaction * const in, KinePhaseSpace_t ps)
{
  double vol = 0;

  const KPhaseSpace & phase_space = in->PhaseSpace();

  switch(ps) {

  case(kPSQ2fE):
  {
    Range1D_t Q2 = phase_space.Limits(kKVQ2);
    vol = Q2.max - Q2.min;
    return vol;
    break;
  }
  case(kPSq2fE):
  {
    Range1D_t q2 = phase_space.Limits(kKVq2);
    vol = q2.min - q2.max;
    return vol;
    break;
  }
  case(kPSWfE):
  {
    Range1D_t W = phase_space.Limits(kKVW);
    vol = W.max - W.min;
    return vol;
    break;
  }
  case(kPSWQ2fE):
  {
    Range1D_t W = phase_space.Limits(kKVW);
    if(W.max<0) return 0;
    const int kNW = 100;
    double dW = (W.max-W.min)/(kNW-1);
    Interaction interaction(*in);
    for(int iw=0; iw<kNW; iw++) {
      interaction.KinePtr()->SetW(W.min + iw*dW);
      Range1D_t Q2 = interaction.PhaseSpace().Q2Lim_W();
      double dQ2 = (Q2.max-Q2.min);
      vol += (dW*dQ2);
    }
    return vol;
    break;
  }
  case(kPSxyfE):
  {
    Range1D_t W = phase_space.Limits(kKVW);
    if(W.max<0) return 0;

    const InitialState & init_state = in->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M();

    const int    kNx = 100;
    const int    kNy = 100;
    const double kminx = controls::kMinX;
    const double kminy = controls::kMinY;
    const double kdx = (controls::kMaxX - kminx) / (kNx-1);
    const double kdy = (controls::kMaxY - kminy) / (kNy-1);
    const double kdV = kdx*kdy;

    double cW=-1, cQ2 = -1;

    Interaction interaction(*in);

    for(int ix=0; ix<kNx; ix++) {
      double x = kminx+ix*kdx;
      for(int iy=0; iy<kNy; iy++) {
         double y = kminy+iy*kdy;

         XYtoWQ2(Ev, M, cW, cQ2, x, y);
         if(!math::IsWithinLimits(cW, W)) continue;

         interaction.KinePtr()->SetW(cW);
         Range1D_t Q2 = interaction.PhaseSpace().Q2Lim_W();
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
  return 0;
}
//____________________________________________________________________________
double genie::utils::kinematics::Jacobian(
 const Interaction * const i, KinePhaseSpace_t fromps, KinePhaseSpace_t tops)
{
// Returns the Jacobian for a kinematical transformation:
// from_ps (x,y,z,...) -> to_ps (u,v,w,...).
//
// Note on the convention used here:
// The differential cross-section d^n sigma / dxdydz..., is transformed to
// d^n sigma / dudvdw... using the Jacobian, computed in this function, as:
// as d^n sigma / dudvdw... = J * d^n sigma / dxdydz...
// where:
//
// dxdxdz... = J * dudvdw...
//
//                    | dx/du   dx/dv   dx/dw  ... |
//                    |                            |
//     d {x,y,z,...}  | dy/du   dy/dv   dy/dw  ... |
// J = ------------ = |                            |
//     d {u,v,w,...}  | dz/du   dz/dv   dz/dw  ... |
//                    |                            |
//                    |  ...     ...    ...    ... |
//
  SLOG("KineLimits", pDEBUG)
       << "Computing Jacobian for transformation: "
       << KinePhaseSpace::AsString(fromps) << " --> "
       << KinePhaseSpace::AsString(tops);

  double J=0;
  bool forward;
  const Kinematics & kine = i->Kine();

  //
  // transformation: {Q2}|E -> {lnQ2}|E
  //
  if ( TransformMatched(fromps,tops,kPSQ2fE,kPSlogQ2fE,forward) )
  {
    J = 1. / kine.Q2();
  }

  //
  // transformation: {QD2}|E -> {Q2}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSQD2fE,kPSQ2fE,forward) )
  {
    J = TMath::Power(1+kine.Q2()/controls::kMQD2,-2)/controls::kMQD2;
  }

  //
  // transformation: {x,y}|E -> {lnx,lny}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSxyfE,kPSlogxlogyfE,forward) )
  {
    J = 1. / (kine.x() * kine.y());
  }

  //
  // transformation: {Q2,y}|E -> {lnQ2,lny}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSQ2yfE,kPSlogQ2logyfE,forward) )
  {
    J = 1. / (kine.Q2() * kine.y());
  }

  //
  // transformation: {Q2,y}|E -> {x,y}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSQ2yfE,kPSxyfE,forward) )
  {
    const InitialState & init_state = i->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M();
    double y  = kine.y();
    J = 2*y*Ev*M;
  }

  //
  // transformation: {W,Q2}|E -> {W,lnQ2}|E
  //
  else if ( TransformMatched(fromps,tops,kPSWQ2fE,kPSWlogQ2fE,forward) )
  {
    J = 1. / kine.Q2();
  }

  //
  // transformation: {W,QD2}|E -> {W,Q2}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSWQD2fE,kPSWQ2fE,forward) )
  {
    J = TMath::Power(1+kine.Q2()/controls::kMQD2,-2)/controls::kMQD2;
  }

  //
  // transformation: {W2,Q2}|E --> {x,y}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSW2Q2fE,kPSxyfE,forward) )
  {
    const InitialState & init_state = i->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M();
    double y  = kine.y();
    J = TMath::Power(2*M*Ev,2) * y;
  }

  //
  // transformation: {W,Q2}|E -> {x,y}|E
  //
  else
  if ( TransformMatched(fromps,tops,kPSWQ2fE,kPSxyfE,forward) )
  {
    const InitialState & init_state = i->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M();
    double y  = kine.y();
    double W  = kine.W();
    J = 2*TMath::Power(M*Ev,2) * y/W;
  }

  // Transformation: {Omegalep,Omegapi}|E -> {Omegalep,Thetapi}|E
  else if ( TransformMatched(fromps,tops,kPSElOlOpifE,kPSElOlTpifE,forward) ) {
    // Use symmetry to turn 4d integral into 3d * 2pi
    J = 2*constants::kPi;
  }

  else {
     SLOG("KineLimits", pFATAL)
       << "*** Can not compute Jacobian for transforming: "
       << KinePhaseSpace::AsString(fromps) << " --> "
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

  // match a->b or b->a
  bool matched = ( (a==inpa&&b==inpb) || (a==inpb&&b==inpa) );

  if(matched) {
    if(a==inpa&&b==inpb) fwd = true;  // matched forward transform: a->b
    else                 fwd = false; // matched reverse transform: b->a
  }
  return matched;
}
//____________________________________________________________________________
// Kinematical Limits:
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelWLim(double Ev, double M, double ml)
{
// Computes W limits for inelastic v interactions
//
  double M2 = TMath::Power(M,2);
  double s  = M2 + 2*M*Ev;
  assert (s>0);

  Range1D_t W;
  W.min  = kNeutronMass + kPhotontest;
  W.max  = TMath::Sqrt(s) - ml;
  if(W.max<=W.min) {
    W.min = -1;
    W.max = -1;
    return W;
  }
  W.min  += controls::kASmallNum;
  W.max  -= controls::kASmallNum;
  return W;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelQ2Lim_W(
     double Ev, double M, double ml, double W, double Q2min_cut)
{
// Computes Q2 limits (>0) @ the input W for inelastic v interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  double M2  = TMath::Power(M,  2.);
  double ml2 = TMath::Power(ml, 2.);
  double W2  = TMath::Power(W,  2.);
  double s   = M2 + 2*M*Ev;

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;
  assert (s>0);

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
  if(Q2.min < Q2min_cut) {Q2.min = Q2min_cut;     }
  if(Q2.max < Q2.min   ) {Q2.min = -1; Q2.max = -1;}

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Inelq2Lim_W(
    double Ev, double M, double ml, double W, double q2min_cut)
{
// Computes q2 (<0) limits @ the input W for inelastic v interactions

  Range1D_t Q2 = utils::kinematics::InelQ2Lim_W(Ev,M,ml,W,-1.*q2min_cut);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelQ2Lim(
    double Ev, double M, double ml, double Q2min_cut)
{
// Computes Q2 (>0) limits irrespective of W for inelastic v interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  Range1D_t W  = utils::kinematics::InelWLim(Ev,M,ml);
  if(W.min<0) return Q2;

  Q2 = utils::kinematics::InelQ2Lim_W(Ev,M,ml,W.min,Q2min_cut);
  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Inelq2Lim(
     double Ev, double M, double ml, double q2min_cut)
{
// Computes Q2 (>0) limits irrespective of W for inelastic v interactions

  Range1D_t Q2 = utils::kinematics::InelQ2Lim(Ev,M,ml,-1.*q2min_cut);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelXLim(double Ev, double M, double ml)

{
// Computes Bjorken x limits for inelastic v interactions

  double M2  = TMath::Power(M, 2.);
  double ml2 = TMath::Power(ml,2.);
  double s   = M2 + 2*M*Ev;

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;
  assert (s>M2);

  Range1D_t x;
  x.min = ml2/(s-M2) + controls::kASmallNum;
  x.max = 1.         - controls::kASmallNum;

  return x;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelYLim(double Ev, double M, double ml)
{
// Computes y limits for inelastic v interactions

  Range1D_t y;
  y.min =  999;
  y.max = -999;

  Range1D_t xl = kinematics::InelXLim(Ev,M,ml);
  assert(xl.min>0 && xl.max>0);

  const unsigned int N=100;
  const double logxmin = TMath::Log10(xl.min);
  const double logxmax = TMath::Log10(xl.max);
  const double dlogx   = (logxmax-logxmin) / (double)(N-1);

  for(unsigned int i=0; i<N; i++) {
    double x = TMath::Power(10, logxmin + i*dlogx);

    Range1D_t y_x = kinematics::InelYLim_X(Ev,M,ml,x);
    if(y_x.min>=0 && y_x.min<=1) y.min = TMath::Min(y.min, y_x.min);
    if(y_x.max>=0 && y_x.max<=1) y.max = TMath::Max(y.max, y_x.max);
  }

  if(y.max >= 0 && y.max <= 1 && y.min >= 0 && y.min <= 1) {
    y.min = TMath::Max(y.min,     controls::kASmallNum);
    y.max = TMath::Min(y.max, 1 - controls::kASmallNum);
  } else {
    y.min = -1;
    y.max = -1;
  }
  SLOG("KineLimits", pDEBUG) << "y  = [" << y.min << ", " << y.max << "]";
  return y;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::InelYLim_X(
                                 double Ev, double M, double ml, double x)

{
// Computes y limits @ the input x for inelastic v interactions

  Range1D_t y;
  y.min = -1;
  y.max = -1;

  double Ev2 = TMath::Power(Ev,2);
  double ml2 = TMath::Power(ml,2);

  SLOG("KineLimits", pDEBUG) << "x  = " << x;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;

  assert (Ev>0);
  assert (x>0&&x<1);

  double a = 0.5 * ml2/(M*Ev*x);
  double b = ml2/Ev2;
  double c = 1 + 0.5*x*M/Ev;
  double d = TMath::Max(0., TMath::Power(1-a,2.) - b);

  double A = 0.5 * (1-a-0.5*b)/c;
  double B = 0.5 * TMath::Sqrt(d)/c;

  y.min = TMath::Max(0., A-B) + controls::kASmallNum;
  y.max = TMath::Min(1., A+B) - controls::kASmallNum;

  return y;
}
//____________________________________________________________________________
// Kinematical Limits for em interactions:
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelWLim(double El, double ml, double M)
{
// Computes W limits for inelastic em interactions
//
  double M2 = TMath::Power(M,2);
  double ml2 = TMath::Power(ml,2); // added lepton mass squared to be used in s calculation
  double s  = M2 + 2*M*El + ml2; // non-negligible mass for em interactions
  assert (s>0);

  Range1D_t W;
  W.min  = kNeutronMass + kPhotontest;
  W.max  = TMath::Sqrt(s) - ml;
  if(W.max<=W.min) {
    W.min = -1;
    W.max = -1;
    return W;
  }
  W.min  += controls::kASmallNum;
  W.max  -= controls::kASmallNum;
  return W;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelQ2Lim_W(
     double El, double ml, double M, double W)
{
// Computes Q2 limits (>0) @ the input W for inelastic em interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  double M2  = TMath::Power(M,  2.);
  double ml2 = TMath::Power(ml, 2.);
  double W2  = TMath::Power(W,  2.);
  double s   = M2 + 2*M*El + ml2; // added lepton mass squared to be used in s calculation (non-negligible mass for em interactions)

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "El = " << El;
  assert (s>0);

  double auxC = 0.5*(s - M2 - ml2)/s; // subtract ml2 to account for the non-negligible mass of the incoming lepton
  double aux1 = s + ml2 - W2;
  double aux2 = aux1*aux1 - 4*s*ml2;

  (aux2 < 0) ? ( aux2 = 0 ) : ( aux2 = TMath::Sqrt(aux2) );

  Q2.max = -ml2 + auxC * (aux1 + aux2); // => 0
  Q2.min = -ml2 + auxC * (aux1 - aux2); // => 0

  // guard against overflows
  Q2.max = TMath::Max(0., Q2.max);
  Q2.min = TMath::Max(0., Q2.min);

  // limit the minimum Q2
  if(Q2.min < utils::kinematics::electromagnetic::kMinQ2Limit) {Q2.min = utils::kinematics::electromagnetic::kMinQ2Limit; } // use the relevant threshold for em scattering 
  if(Q2.max < Q2.min   ) {Q2.min = -1; Q2.max = -1;}

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::Inelq2Lim_W(
    double El, double ml, double M, double W)
{
// Computes q2 (<0) limits @ the input W for inelastic em interactions

  Range1D_t Q2 = utils::kinematics::electromagnetic::InelQ2Lim_W(El,ml,M,W);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelQ2Lim(
    double El, double ml, double M)
{
// Computes Q2 (>0) limits irrespective of W for inelastic em interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  Range1D_t W  = utils::kinematics::electromagnetic::InelWLim(El,ml,M);
  if(W.min<0) return Q2;

  Q2 = utils::kinematics::electromagnetic::InelQ2Lim_W(El,ml,M,W.min);
  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::Inelq2Lim(
     double El, double ml, double M)
{
// Computes Q2 (>0) limits irrespective of W for inelastic em interactions

  Range1D_t Q2 = utils::kinematics::electromagnetic::InelQ2Lim(El,ml,M);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelXLim(double El, double ml, double M)

{
// Computes Bjorken x limits for inelastic em interactions

  double M2  = TMath::Power(M, 2.);
  double ml2 = TMath::Power(ml,2.);
  double s   = M2 + 2*M*El + ml2; // added lepton mass squared to be used in s calculation (non-negligible mass for em interactions)

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "El = " << El;
  assert (s > M2 + ml2); // added lepton mass squared to be used in s calculation (non-negligible mass for em interactions)

  Range1D_t x;
  x.min = ml2/(s - M2 - ml2) + controls::kASmallNum; // subtracted lepton mass squared to be used in s calculation (non-negligible mass for em interactions)
  x.max = 1. - controls::kASmallNum;

  return x;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelYLim(double El, double ml, double M)
{
// Computes y limits for inelastic v interactions

  Range1D_t y;
  y.min =  999;
  y.max = -999;

  Range1D_t xl = kinematics::electromagnetic::InelXLim(El,ml,M);
  assert(xl.min>0 && xl.max>0);

  const unsigned int N=100;
  const double logxmin = TMath::Log10(xl.min);
  const double logxmax = TMath::Log10(xl.max);
  const double dlogx   = (logxmax-logxmin) / (double)(N-1);

  for(unsigned int i=0; i<N; i++) {
    double x = TMath::Power(10, logxmin + i*dlogx);

    Range1D_t y_x = kinematics::electromagnetic::InelYLim_X(El,ml,M,x);
    if(y_x.min>=0 && y_x.min<=1) y.min = TMath::Min(y.min, y_x.min);
    if(y_x.max>=0 && y_x.max<=1) y.max = TMath::Max(y.max, y_x.max);
  }

  if(y.max >= 0 && y.max <= 1 && y.min >= 0 && y.min <= 1) {
    y.min = TMath::Max(y.min,     controls::kASmallNum);
    y.max = TMath::Min(y.max, 1 - controls::kASmallNum);
  } else {
    y.min = -1;
    y.max = -1;
  }
  SLOG("KineLimits", pDEBUG) << "y  = [" << y.min << ", " << y.max << "]";
  return y;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::electromagnetic::InelYLim_X(
                                 double El, double ml, double M, double x)

{
// Computes y limits @ the input x for inelastic em interactions

  Range1D_t y;
  y.min = -1;
  y.max = -1;

  double El2 = TMath::Power(El,2);
  double ml2 = TMath::Power(ml,2);

  SLOG("KineLimits", pDEBUG) << "x  = " << x;
  SLOG("KineLimits", pDEBUG) << "El = " << El;

  assert (El>0);
  assert (x>0&&x<1);

  double a = 0.5 * ml2/(M*El*x);
  double b = ml2/El2;
  double c = 1 + 0.5*x*M/El;
  double d = TMath::Max(0., TMath::Power(1-a,2.) - b);

  double A = 0.5 * (1-a-0.5*b)/c;
  double B = 0.5 * TMath::Sqrt(d)/c;

  y.min = TMath::Max(0., A-B) + controls::kASmallNum;
  y.max = TMath::Min(1., A+B) - controls::kASmallNum;

  return y;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohXLim(void)
{
// Computes x limits for coherent v interactions

  Range1D_t x(controls::kASmallNum, 1.-controls::kASmallNum);
  return x;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohQ2Lim(double Mn, double mpi, double mlep, double Ev)
{
  // The expressions for Q^2 min appears in PRD 74, 054007 (2006) by
  // Kartavtsev, Paschos, and Gounaris

  Range1D_t Q2;
  Q2.min = 0.0;
  Q2.max = std::numeric_limits<double>::max();  // Value must be overriden in user options

  double Mn2 = Mn * Mn;
  double mlep2 = mlep * mlep;
  double s = Mn2 + 2.0 * Mn * Ev;
  double W2min = CohW2Min(Mn, mpi);

  // Looks like Q2min = A * B - C, where A, B, and C are complicated
  double a = 1.0;
  double b = mlep2 / s;
  double c = W2min / s;
  double lambda = a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c;
  if (lambda > 0) {
    double A = (s - Mn * Mn) / 2.0;
    double B = 1 - TMath::Sqrt(lambda);
    double C = 0.5 * (W2min + mlep2 - Mn2 * (W2min - mlep2) / s );
    if (A * B - C < 0) {
      SLOG("KineLimits", pERROR)
        << "Q2 kinematic limits calculation failed for CohQ2Lim. "
        << "Assuming Q2min = 0.0";
    }
    Q2.min = TMath::Max(0., A * B - C);
  } else {
    SLOG("KineLimits", pERROR)
      << "Q2 kinematic limits calculation failed for CohQ2Lim. "
      << "Assuming Q2min = 0.0";
  }

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Cohq2Lim(double Mn, double mpi, double mlep, double Ev)
{
  Range1D_t Q2 = utils::kinematics::CohQ2Lim(Mn, mpi, mlep, Ev);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohW2Lim(double Mn, double mpi, double mlep,
    double Ev, double Q2)
{
  // These expressions for W^2 min and max appear in PRD 74, 054007 (2006) by
  // Kartavtsev, Paschos, and Gounaris

  Range1D_t W2l;
  W2l.min = -1;
  W2l.max = -1;

  double s = Mn * Mn + 2.0 * Mn * Ev;
  double Mnterm = 1 - Mn * Mn / s;
  double Mlterm = 1 - mlep * mlep / s;
  // Here T1, T2 are generically "term 1" and "term 2" in a long expression
  double T1 = 0.25 * s * s * Mnterm * Mnterm * Mlterm;
  double T2 = Q2 - (0.5 * s * Mnterm) + (0.5 * mlep * mlep * Mnterm);

  W2l.min = CohW2Min(Mn, mpi);
  W2l.max = (T1 - T2 * T2 ) *
    (1.0 / Mnterm) *
    (1.0 / (Q2 + mlep * mlep));

  return W2l;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohNuLim(double W2min, double W2max,
    double Q2, double Mn, double xsi)
{
  Range1D_t nul;
  nul.min = -1;
  nul.max = -1;

  double nu_min = (W2min + Q2 - Mn * Mn) / (2.0 * Mn);
  double nu_max = (W2max + Q2 - Mn * Mn) / (2.0 * Mn);
  double xsiQ = xsi * TMath::Sqrt(Q2);

  nul.min = (xsiQ > nu_min) ? xsiQ : nu_min;
  nul.max = nu_max;

  return nul;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohYLim(double Mn, double mpi, double mlep,
    double Ev, double Q2, double xsi)
{
  Range1D_t ylim;
  ylim.min = -1;
  ylim.max = -1;

  Range1D_t W2lim = genie::utils::kinematics::CohW2Lim(Mn, mpi, mlep, Ev, Q2);
  if (W2lim.min > W2lim.max) {
    LOG("KineLimits", pDEBUG)
      << "Kinematically forbidden region in CohYLim. W2min = " << W2lim.min
      << "; W2max =" << W2lim.max;
    LOG("KineLimits", pDEBUG)
      << "  Mn = " << Mn << "; mpi = " << mpi << "; mlep = "
      << mlep << "; Ev = " << Ev << "; Q2 = " << Q2;
    return ylim;
  }
  Range1D_t nulim = genie::utils::kinematics::CohNuLim(W2lim.min, W2lim.max,
      Q2, Mn, xsi);
  ylim.min = nulim.min / Ev;
  ylim.max = nulim.max / Ev;

  return ylim;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::CohYLim(double EvL, double ml)
{
// Computes y limits for coherent v interactions

  Range1D_t y(kPionMass/EvL + controls::kASmallNum,
              1.-ml/EvL - controls::kASmallNum);
  return y;
}
//____________________________________________________________________________
// Helpers for kinematic limits
//____________________________________________________________________________
double genie::utils::kinematics::CohW2Min(double Mn, double mpi)
{
  // These expressions for W^2 min and max appear in PRD 74, 054007 (2006) by
  // Kartavtsev, Paschos, and Gounaris

  return (Mn + mpi) * (Mn + mpi);
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkWLim(double Ev, double M, double ml)
{
// Computes W limits for inelastic v interactions
//
  double M2 = TMath::Power(M,2);
  double ml2 = TMath::Power(ml,2);
  double s  = M2 + 2*M*Ev + ml2;
  assert (s>0);

  Range1D_t W;
  W.min  = kNeutronMass + kPhotontest;
//  W.min  = kNeutronMass + kPionMass;
  W.max  = TMath::Sqrt(s) - ml;
  if(W.max<=W.min) {
    W.min = -1;
    W.max = -1;
    return W;
  }
  W.min  += controls::kASmallNum;
  W.max  -= controls::kASmallNum;
  return W;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkQ2Lim_W(
     double Ev, double M, double ml, double W, double Q2min_cut)
{
// Computes Q2 limits (>0) @ the input W for inelastic v interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  double M2  = TMath::Power(M,  2.);
  double ml2 = TMath::Power(ml, 2.);
  double W2  = TMath::Power(W,  2.);
  double s   = M2 + 2*M*Ev + ml2;
  assert(s > 0);
  double sqs = TMath::Sqrt(s);
  double E1CM = (s + ml2 - M2) / (2.*sqs);
  double p1CM = TMath::Max(0., E1CM*E1CM - ml2);
  p1CM = TMath::Sqrt(p1CM);
  double E3CM = (s + ml2 - W2) / (2.*sqs);
  double p3CM = TMath::Max(0., E3CM*E3CM - ml2);
  p3CM = TMath::Sqrt(p3CM);

  SLOG("KineLimits", pDEBUG) << "s  = " << s;
  SLOG("KineLimits", pDEBUG) << "Ev = " << Ev;
  SLOG("KineLimits", pDEBUG) << "M = " << M;
  SLOG("KineLimits", pDEBUG) << "W = " << W;
  SLOG("KineLimits", pDEBUG) << "E1_CM = " << E1CM;
  SLOG("KineLimits", pDEBUG) << "p1_CM = " << p1CM;
  SLOG("KineLimits", pDEBUG) << "E3_CM = " << E3CM;
  SLOG("KineLimits", pDEBUG) << "p3_CM = " << p3CM;

  Q2.min = TMath::Power(p3CM - p1CM,2) - TMath::Power((W2 - M2) / (2.*sqs),2);
  Q2.max = TMath::Power(p3CM + p1CM,2) - TMath::Power((W2 - M2) / (2.*sqs),2);

  SLOG("KineLimits", pDEBUG) << "Nominal Q^2 limits: " << Q2.min << " , " << Q2.max;
  // guard against overflows
  Q2.max = TMath::Max(0., Q2.max);
  Q2.min = TMath::Max(0., Q2.min);

  // limit the minimum Q2
  if(Q2.min < Q2min_cut) {Q2.min = Q2min_cut;     }
  if(Q2.max < Q2.min   ) {Q2.min = -1; Q2.max = -1;}

  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Darkq2Lim_W(
    double Ev, double M, double ml, double W, double q2min_cut)
{
// Computes q2 (<0) limits @ the input W for inelastic v interactions

  Range1D_t Q2 = utils::kinematics::DarkQ2Lim_W(Ev,M,ml,W,-1.*q2min_cut);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkQ2Lim(
    double Ev, double M, double ml, double Q2min_cut)
{
// Computes Q2 (>0) limits irrespective of W for inelastic v interactions

  Range1D_t Q2;
  Q2.min = -1;
  Q2.max = -1;

  Range1D_t W  = utils::kinematics::DarkWLim(Ev,M,ml);
  if(W.min<0) return Q2;

  Q2 = utils::kinematics::DarkQ2Lim_W(Ev,M,ml,W.min,Q2min_cut);
  return Q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::Darkq2Lim(
     double Ev, double M, double ml, double q2min_cut)
{
// Computes Q2 (>0) limits irrespective of W for inelastic v interactions

  Range1D_t Q2 = utils::kinematics::DarkQ2Lim(Ev,M,ml,-1.*q2min_cut);
  Range1D_t q2;
  q2.min = - Q2.max;
  q2.max = - Q2.min;
  return q2;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkXLim(double Ev, double M, double ml)

{
// Computes Bjorken x limits for inelastic interactions
// For the dark matter case, it is relatively straightforward
  Range1D_t Wl = utils::kinematics::DarkWLim(Ev, M, ml);
  double Wmin = Wl.min;
  double W2min = Wmin*Wmin;
  SLOG("KineLimits", pDEBUG) << "W^2_min = " << W2min;
  Range1D_t Q2l = utils::kinematics::DarkQ2Lim_W(Ev, M, ml, Wmin);
  SLOG("KineLimits", pDEBUG) << "Q^2 range : " << Q2l.min << " , " << Q2l.max;
  double M2 = M*M;
  Range1D_t x;
  x.min = Q2l.min / (Q2l.min + W2min - M2);
  x.max = Q2l.max / (Q2l.max + W2min - M2);

  SLOG("KineLimits", pDEBUG) << "x  = [" << x.min << ", " << x.max << "]";
  return x;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkYLim(double Ev, double M, double ml)
{
  // For dark inelastic scattering, can compute exactly and is much simpler
  Range1D_t Wl = utils::kinematics::DarkWLim(Ev, M, ml);
  double Wmin = Wl.min;
  double W2min = Wmin*Wmin;
  Range1D_t Q2l = utils::kinematics::DarkQ2Lim_W(Ev, M, ml, Wmin);
  double M2 = M*M;
  Range1D_t y;
  y.min = (Q2l.min + W2min - M2) / (2*Ev*M);
  y.max = (Q2l.max + W2min - M2) / (2*Ev*M);

  SLOG("KineLimits", pDEBUG) << "y  = [" << y.min << ", " << y.max << "]";
  return y;
}
//____________________________________________________________________________
Range1D_t genie::utils::kinematics::DarkYLim_X(
                                 double Ev, double M, double ml, double x)

{
  // Computes y limits @ the input x for inelastic interactions
  // We hit y_min when W = W_min and y_max when Q^2 = Q_min^2 or Q^2 = Q_max^2

  Range1D_t y;
  y.min = -1;
  y.max = -1;

  Range1D_t Wl = utils::kinematics::DarkWLim(Ev, M, ml);
  double Wmin = Wl.min;
  double W2min = Wmin*Wmin;
  double M2 = M*M;
  double ml2 = ml*ml;
  y.min = (W2min - M2) / (1.-x) / (2*Ev*M);
  y.max = 2.* M * x *(Ev*Ev - ml2) / Ev / (2. * M * Ev * x + M2 * x * x + ml2);

  return y;
}
//____________________________________________________________________________
// Kinematical Transformations:
//____________________________________________________________________________
double genie::utils::kinematics::Q2toQD2(double Q2)
{
// Q2 -> QD2 transformation. QD2 takes out the dipole form of the form factors
// making the differential cross section to be flatter and speeding up the
// kinematical selection.

  assert(Q2>0);
  return TMath::Power(1+Q2/controls::kMQD2, -1);
}
//____________________________________________________________________________
double genie::utils::kinematics::QD2toQ2(double QD2)
{
  assert(QD2>0);
  return controls::kMQD2*(1/QD2-1);
}
//____________________________________________________________________________
double genie::utils::kinematics::Q2(const Interaction * const interaction)
{
// Get Q^2 from kinematics object
// If Q^2 is not set but x,y are, then compute Q2 from x,y

  const Kinematics & kinematics = interaction->Kine();

  if (kinematics.KVSet(kKVQ2) || kinematics.KVSet(kKVq2)) {
    double Q2 = kinematics.Q2();
    return Q2;
  }
  if (kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->InitState();
    double Mn = init_state.Tgt().HitNucP4Ptr()->M(); // can be off m/shell
    double x  = kinematics.x();
    double y  = kinematics.y();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double Q2 = 2*Mn*Ev*x*y;
    return Q2;
  }
  SLOG("KineLimits", pERROR) << "Couldn't compute Q^2 for \n"<< *interaction;
  return 0;
}
//____________________________________________________________________________
double genie::utils::kinematics::W(const Interaction * const interaction)
{
  const ProcessInfo & process_info = interaction->ProcInfo();

  if(process_info.IsQuasiElastic()) {
    // hadronic inv. mass is equal to the recoil nucleon on-shell mass
    int rpdgc = interaction->RecoilNucleonPdg();
    double M = PDGLibrary::Instance()->Find(rpdgc)->Mass();
    return M;
  }

  const Kinematics & kinematics = interaction->Kine();
  if(kinematics.KVSet(kKVW)) {
    double W = kinematics.W();
    return W;
  }
  if(kinematics.KVSet(kKVx) && kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M();
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
//___________________________________________________________________________
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
double genie::utils::kinematics::Q2YtoX(
                                    double Ev, double M, double Q2, double y)
{
// Converts (Q^2,y) => x
// Ev is the neutrino energy at the struck nucleon rest frame
// M is the nucleon mass - it does not need to be on the mass shell
  assert(Ev > 0. && M  > 0. && Q2 > 0. && y  > 0.);

  double x = Q2 / (2. * y * M * Ev);

  LOG("KineLimits", pDEBUG) << "(Ev=" << Ev << ",Q2=" << Q2
    << ",y=" << y << ",M=" << M << ") => x=" << x;

  return x;
}
//___________________________________________________________________________
bool genie::utils::kinematics::IsAboveCharmThreshold(
                                   double x, double Q2, double M, double mc)
{
// x : scaling variable (plain or modified)
// Q2: momentum transfer
// M : hit nucleon "mass" (nucleon can be off the mass shell)
// mc: charm mass
//
  double M2   = TMath::Power(M,2);
  double v    = 0.5*Q2/(M*x);
  double W2   = TMath::Max(0., M2+2*M*v-Q2);
  double W    = TMath::Sqrt(W2);
  double Wmin = M + kLightestChmHad;
  double xc   = utils::kinematics::SlowRescalingVar(x,Q2,M,mc);

  if(xc>=1 || W<=Wmin) return false;
  else                 return true;
}
//____________________________________________________________________________
double genie::utils::kinematics::SlowRescalingVar(
                                   double x, double Q2, double M, double mc)
{
// x : scaling variable (plain or modified)
// Q2: momentum transfer
// M : hit nucleon "mass" (nucleon can be off the mass shell)
// mc: charm mass

  double mc2 = TMath::Power(mc,2);
  double v   = 0.5*Q2/(M*x);
  double xc  = x + 0.5*mc2/(M*v);
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
//___________________________________________________________________________
void genie::utils::kinematics::UpdateWQ2FromXY(const Interaction * in)
{
  Kinematics * kine = in->KinePtr();

  if(kine->KVSet(kKVx) && kine->KVSet(kKVy)) {
    const InitialState & init_state = in->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M(); // can be off mass shell
    double x  = kine->x();
    double y  = kine->y();

    double Q2=-1,W=-1;
    kinematics::XYtoWQ2(Ev,M,W,Q2,x,y);
    kine->SetQ2(Q2);
    kine->SetW(W);
  }
}
//___________________________________________________________________________
void genie::utils::kinematics::UpdateXYFromWQ2(const Interaction * in)
{
  Kinematics * kine = in->KinePtr();

  if(kine->KVSet(kKVW) && kine->KVSet(kKVQ2)) {
    const InitialState & init_state = in->InitState();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double M  = init_state.Tgt().HitNucP4Ptr()->M(); // can be off mass shell
    double W  = kine->W();
    double Q2 = kine->Q2();

    double x=-1,y=-1;
    kinematics::WQ2toXY(Ev,M,W,Q2,x,y);
    kine->Setx(x);
    kine->Sety(y);
  }
}
//___________________________________________________________________________
void genie::utils::kinematics::UpdateXFromQ2Y(const Interaction * in)
{
  Kinematics * kine = in->KinePtr();

  if(kine->KVSet(kKVy) && kine->KVSet(kKVQ2)) {

    const ProcessInfo &  pi = in->ProcInfo();
    const InitialState & init_state = in->InitState();
    double M = 0.0;
    double Ev = 0.0;

    if (pi.IsCoherent()) {
      M = in->InitState().Tgt().Mass(); // nucleus mass
      Ev = init_state.ProbeE(kRfLab);
    }
    else {
      M = in->InitState().Tgt().HitNucP4Ptr()->M(); //can be off m/shell
      Ev = init_state.ProbeE(kRfHitNucRest);
    }

    double y  = kine->y();
    double Q2 = kine->Q2();

    double x = kinematics::Q2YtoX(Ev,M,Q2,y);
    kine->Setx(x);
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
  //double E     = par[4]; // neutrino energy

//  return 3*xsmax; ////////////////NOTE

  xsmax*=5;

  double func = 0;

  if(Wmax > mD) {
     // -- if the resonance mass is within the kinematical range then
     // -- the envelope is defined as a plateau above the resonance peak,
     // -- a steeply falling leading edge (low-W side) and a more slowly
     // -- falling trailing edge (high-W side)

     double hwfe = mD+2*gD; // high W falling edge
     double lwfe = mD-gD/2; // low  W falling edge

     if(W < lwfe) {
       //low-W falling edge
       func = xsmax / (1 + 5* TMath::Power((W-lwfe)/gD,2));
     } else if (W > hwfe) {
       //high-W falling edge
       func = xsmax / (1 + TMath::Power((W-hwfe)/gD,2));
     } else {
       // plateau
       func = xsmax;
     }
  } else {
     // -- if the resonance mass is above the kinematical range then the
     // -- envelope is a small plateau just bellow Wmax and a falling edge
     // -- at lower W

     double plateau_edge = Wmax-0.1;

     if (W > plateau_edge) {
       // plateau
       func = xsmax;
     } else {
       //low-W falling edge
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
  double yb = x[1];   // y

  //-- parameters
  double xsmax = 3*par[0]; // safety factor * max cross section in (x,y)
  double Ev    = par[1]; // neutrino energy;

  if(yb<0.|| yb>1.) return 0.;
  if(xb<0.|| xb>1.) return 0.;

  if(Ev<1) return xsmax;
  if(xb/Ev<1E-4 && yb>0.95) return 5*xsmax;

  double func = 0;
  double xp   = 0.1;
  double yp   = (Ev>2.5) ? 2.5/Ev : 1;

  if(xb>xp) {
    double xs0=0;
    if(yb<yp) {
      xs0 = xsmax;
    } else {
      xs0 = xsmax - (yb-yp)*xsmax;
    }
    double d = TMath::Power( (xb-xp)/0.075, 2);
    func = xs0/(1 + d);
  } else {
    if(yb>yp) {
       func = xsmax - (yb-yp)*xsmax;
    } else {
       func = xsmax;
    }
  }
  return func;
}
//___________________________________________________________________________
