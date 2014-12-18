//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Removed the hardcoded Q2min cut to accomododate the group extending GENIE
   to ~MeV range energies (They require the Q2min cut to go down from 1E-4
   GeV^2 to 1E-10 GeV^2). Now all methods computing Q2 limits accept the cut
   as argument. The argument has a default value (controls::kMinQ2Limit) that 
   leaves the intermediate energy (~GeV) simulations unaffected.
   Removed 'using namespace genie::controls' and added the controls namespace
   qualifier at all relevent consts to clarify their source.
   
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"

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
#ifdef __GENIE_VHE_ENABLED__ 
    const double kminx = controls::kMinX_VHE;	
    const double kminy = controls::kMinY_VHE;	
#else
    const double kminx = controls::kMinX;	
    const double kminy = controls::kMinY;	
#endif
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
  W.min  = kNeutronMass + kPionMass;
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
Range1D_t genie::utils::kinematics::CohXLim(void)
{
// Computes x limits for coherent v interactions

  Range1D_t x(controls::kASmallNum, 1.-controls::kASmallNum);
  return x;
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


