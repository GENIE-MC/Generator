//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HELepton/XSection/Born.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Conventions/Constants.h"

#include <TMath.h>

#include <iostream>

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
Born::Born() 
{

  fGw      = PDGLibrary::Instance()->Find(kPdgWM)->Width();
  fGz      = PDGLibrary::Instance()->Find(kPdgZ0)->Width();
  fmw2c    = TComplex(kMw2,fGw*kMw);
  fmz2c    = TComplex(kMz2,fGz*kMz);
  TComplex rat = fmw2c/fmz2c;
  fsw2     = TComplex(1.-rat.Re(),-rat.Im());
  fcw2     = 1.-fsw2;
  falpha   = TMath::Sqrt(2.)*kGF/kPi * fmw2c * fsw2;
  
  fgLe     = (-1./2.+fsw2)*TMath::Sqrt(1./(fsw2*fcw2) );
  fgRe     = TMath::Sqrt(fsw2/fcw2);
  fgLnu    = 1./2./TMath::Sqrt(fsw2/fcw2);

}
//____________________________________________________________________________
Born::~Born()
{

}
//____________________________________________________________________________
double Born::PXSecCCR(double s, double t, double mlin2, double mlout2)
{
/*
nu \  W.  / nu
    ------
 l /      \ l
*/

  TComplex prop = falpha/fsw2/(s-fmw2c);

  return (t-mlout2)*(t-mlin2) * prop.Rho2();

}
//____________________________________________________________________________
double Born::PXSecCCV(double s, double t, double mlin2, double mlout2)
{
/*
nu -------- l
      | W
 l -------- nu
*/

  return 0.;

}
//____________________________________________________________________________
double Born::PXSecNCV(double s, double t, double mlin2, double mlout2)
{
/*
nu -------- nu
      | Z
 l -------- l
*/

  double u = mlin2 + mlout2 - s - t;
  
  TComplex a = 4.*fgLnu*fgRe/(u-fmz2c);
  TComplex b = 2.*fgLnu*fgLe/(u-fmz2c);   

  return falpha.Rho2() * ( (s-mlout2)*(s-mlin2)*a.Rho2() + (t-mlout2)*(t-mlin2)*b.Rho2() );

}
//____________________________________________________________________________
double Born::PXSecCCRNC(double s, double t, double mlin2, double mlout2)
{
/*
nu \  W.  / nu     nu -------- nu
    ------      +        | Z
 l /      \ l       l -------- l
*/

  double u = mlin2 + mlout2 - s - t;
  
  TComplex a = 4.*fgLnu*fgRe/(u-fmz2c);
  TComplex b = 2.*fgLnu*fgLe/(u-fmz2c)+1./fsw2/(s-fmw2c);   

  return falpha.Rho2() * ( (s-mlout2)*(s-mlin2)*a.Rho2() + (t-mlout2)*(t-mlin2)*b.Rho2() );

}
//____________________________________________________________________________
double Born::PXSecCCVNC(double s, double t, double mlin2, double mlout2)
{
/*
nu -------- l     nu -------- nu
      | W       +        | Z
 l -------- nu       l -------- l
*/

  return 0.;

}
//____________________________________________________________________________
double Born::PXSecPhoton(double s, double t, double mlout2)
{

  double u = kMw2 + mlout2 - s - t;
  
  double ME = 8*kPi2*falpha.Rho2()*s/kMw2/fsw2.Re()/TMath::Power(mlout2 - t,2)/TMath::Power(kMw2 - u,2)* 
        ( -2*(kMw2 - u)*TMath::Power(mlout2,3) - 2*TMath::Power(mlout2,2)*(-2*kMw2*u + u*(s + u) + TMath::Power(kMw2,2)) 
          + mlout2*(-(kMw2*u*(4*s + 5*u)) + (s + u)*TMath::Power(kMw2,2) + 3*TMath::Power(kMw2,3) + (s + u)*TMath::Power(u,2)) 
          + 2*kMw2*((3*s + u)*TMath::Power(kMw2,2) - TMath::Power(kMw2,3) + 4*u*TMath::Power(s,2) + 2*TMath::Power(s,3) + 3*s*TMath::Power(u,2) 
          + TMath::Power(u,3) - kMw2*TMath::Power(2*s + u,2)) );
  
  return TMath::Max(0.,ME);

}
//____________________________________________________________________________
double Born::PXSecPhoton_T(double s12, double s13, double Q2, double ml2){
  double ME2 = 0.0;
  ME2 = (4*falpha.Rho2()*kPi2*(TMath::Power(ml2,4)*s12*(2*TMath::Power(kMw2,2)*TMath::Power(s12,2) - 2*kMw2*Q2*s12*s13 + TMath::Power(Q2,2)*s13*(-s12 + s13)) + 
     TMath::Power(ml2,3)*(-2*TMath::Power(kMw2,3)*TMath::Power(s12,3) + TMath::Power(Q2,2)*(s12 - s13)*s13*(-(Q2*s12) + TMath::Power(s12,2) + Q2*s13 - 3*s12*s13) + 
        2*TMath::Power(kMw2,2)*TMath::Power(s12,2)*(3*Q2*s12 - 2*TMath::Power(s12,2) + 2*Q2*s13 + s12*s13) + 
        kMw2*Q2*s12*(Q2*TMath::Power(s12,2) - Q2*TMath::Power(s13,2) - 2*s12*TMath::Power(s13,2))) + 
     TMath::Power(ml2,2)*(-6*TMath::Power(kMw2,4)*TMath::Power(s12,3) + 2*kMw2*Q2*
         (-(Q2*s12*(s12 - 3*s13)) + TMath::Power(Q2,2)*(s12 - s13) + TMath::Power(s12,2)*(s12 - s13))*(s12 - s13)*s13 + 
        TMath::Power(Q2,2)*(s12 - s13)*TMath::Power(s13,2)*(-2*Q2*s12 + 2*TMath::Power(s12,2) + 2*Q2*s13 - 3*s12*s13) + 
        2*TMath::Power(kMw2,3)*TMath::Power(s12,2)*(2*TMath::Power(s12,2) - 3*Q2*(2*s12 + s13)) + 
        TMath::Power(kMw2,2)*s12*(2*TMath::Power(s12,2)*TMath::Power(s12 - s13,2) + TMath::Power(Q2,2)*(s12 - s13)*s13 - 
           2*Q2*s12*(TMath::Power(s12,2) - 8*s12*s13 - 2*TMath::Power(s13,2)))) - 
     2*TMath::Power(kMw2,2)*TMath::Power(s12,2)*(2*TMath::Power(kMw2,4)*s12 + TMath::Power(Q2,2)*TMath::Power(s13,2)*(-s12 + s13) - 
        2*TMath::Power(kMw2,3)*(2*TMath::Power(s12,2) - Q2*s13 + s12*s13) + 
        TMath::Power(kMw2,2)*(4*Q2*s12*(s12 - 2*s13) + TMath::Power(Q2,2)*(-s12 + s13) + 2*s12*TMath::Power(s12 + s13,2)) + 
        2*kMw2*s13*(TMath::Power(Q2,2)*(s12 - s13) - s12*(TMath::Power(s12,2) + TMath::Power(s13,2)) + Q2*(-TMath::Power(s12,2) + 2*s12*s13 + TMath::Power(s13,2)))) + 
     ml2*(10*TMath::Power(kMw2,5)*TMath::Power(s12,3) - TMath::Power(Q2,2)*(Q2 - s12)*TMath::Power(s12 - s13,2)*TMath::Power(s13,3) + 
        2*TMath::Power(kMw2,4)*TMath::Power(s12,2)*(3*Q2*s12 - 4*TMath::Power(s12,2) + 4*Q2*s13 - 3*s12*s13) + 
        kMw2*Q2*(s12 - s13)*TMath::Power(s13,2)*(2*TMath::Power(Q2,2)*(s12 - s13) + 2*TMath::Power(s12,2)*(s12 - s13) + Q2*s12*(-3*s12 + 5*s13)) + 
        TMath::Power(kMw2,3)*s12*(2*Q2*s12*(5*TMath::Power(s12,2) - 16*s12*s13 - TMath::Power(s13,2)) + 
           TMath::Power(Q2,2)*(-3*TMath::Power(s12,2) + 2*s12*s13 + TMath::Power(s13,2)) + 2*TMath::Power(s12,2)*(TMath::Power(s12,2) + 6*s12*s13 + TMath::Power(s13,2))) - 
        TMath::Power(kMw2,2)*s13*(TMath::Power(Q2,3)*TMath::Power(s12 - s13,2) - 2*TMath::Power(s12,3)*TMath::Power(s12 - s13,2) + 
           TMath::Power(Q2,2)*s12*(-5*TMath::Power(s12,2) + 8*s12*s13 - 3*TMath::Power(s13,2)) + 
           2*Q2*TMath::Power(s12,2)*(3*TMath::Power(s12,2) - 9*s12*s13 + 2*TMath::Power(s13,2))))))/(TMath::Power(kMw2,3)*TMath::Power(s12,2)*TMath::Power(ml2 - kMw2 + s13,2)*TMath::Power(ml2 - kMw2 - s12 + s13,2)*fsw2.Re());
  return TMath::Max(0.,ME2);
}
//____________________________________________________________________________
double Born::PXSecPhoton_L(double s12, double s13, double Q2, double ml2){
  double ME2 = 0.0;
  ME2 = 2*falpha.Rho2()*Q2*TMath::Power(kMw2,-3)*kPi2*TMath::Power(s12,-2)*TMath::Power(ml2 - kMw2 + s13,-2)*TMath::Power(ml2 - kMw2 - s12 + s13,-2)*
 ((s12 - s13)*TMath::Power(ml2,5)*TMath::Power(s12,2) - 2*s12*TMath::Power(ml2,4)*(2*kMw2*TMath::Power(s12,2) - (Q2 - s12)*(-3*s12*s13 + TMath::Power(s12,2) + 2*TMath::Power(s13,2))) + 
   TMath::Power(ml2,3)*(-2*(s12 - s13)*TMath::Power(kMw2,2)*TMath::Power(s12,2) + 
      2*kMw2*s12*(Q2*(9*s12*s13 - 5*TMath::Power(s12,2) - 2*TMath::Power(s13,2)) + s12*(-11*s12*s13 + 3*TMath::Power(s12,2) + 4*TMath::Power(s13,2))) + 
      (s12 - s13)*(TMath::Power(Q2,2)*TMath::Power(s12 - 2*s13,2) + TMath::Power(s12,2)*(-6*s12*s13 + TMath::Power(s12,2) + 6*TMath::Power(s13,2)) - 
         2*Q2*s12*(-5*s12*s13 + TMath::Power(s12,2) + 6*TMath::Power(s13,2)))) + 
   2*TMath::Power(ml2,2)*(4*(2*s12 - s13)*TMath::Power(kMw2,3)*TMath::Power(s12,2) + 
      (Q2 - s12)*s13*(Q2*(s12 - 2*s13) + s12*(-s12 + s13))*(-3*s12*s13 + TMath::Power(s12,2) + 2*TMath::Power(s13,2)) + 
      s12*TMath::Power(kMw2,2)*(Q2*(-11*s12*s13 + 7*TMath::Power(s12,2) - 2*TMath::Power(s13,2)) + s12*(15*s12*s13 - 11*TMath::Power(s12,2) + 2*TMath::Power(s13,2))) - 
      kMw2*((s12 - s13)*TMath::Power(Q2,2)*TMath::Power(s12 - 2*s13,2) + TMath::Power(s12,2)*(-11*s13*TMath::Power(s12,2) + TMath::Power(s12,3) + 20*s12*TMath::Power(s13,2) - 8*TMath::Power(s13,3)) - 
         2*Q2*s12*(-10*s13*TMath::Power(s12,2) + TMath::Power(s12,3) + 17*s12*TMath::Power(s13,2) - 6*TMath::Power(s13,3)))) + 
   4*TMath::Power(kMw2,2)*TMath::Power(s12,2)*(-(s13*(Q2 - 7*s12 + 4*s13)*TMath::Power(kMw2,2)) + (s12 - 2*s13)*TMath::Power(kMw2,3) + kMw2*(2*Q2 - s12 - 2*s13)*TMath::Power(s13,2) + 
      (-Q2 + s12)*TMath::Power(s13,3)) + ml2*(-15*(s12 - s13)*TMath::Power(kMw2,4)*TMath::Power(s12,2) + 
      2*s12*TMath::Power(kMw2,3)*(Q2*(7*s12*s13 - 3*TMath::Power(s12,2) + 2*TMath::Power(s13,2)) + s12*(-21*s12*s13 + 9*TMath::Power(s12,2) + 4*TMath::Power(s13,2))) + 
      TMath::Power(kMw2,2)*((s12 - s13)*TMath::Power(Q2,2)*TMath::Power(s12 - 2*s13,2) + 
         TMath::Power(s12,2)*(-31*s13*TMath::Power(s12,2) + TMath::Power(s12,3) + 60*s12*TMath::Power(s13,2) - 14*TMath::Power(s13,3)) - 
         2*Q2*s12*(-14*s13*TMath::Power(s12,2) + TMath::Power(s12,3) + 27*s12*TMath::Power(s13,2) - 6*TMath::Power(s13,3))) - 
      2*kMw2*s13*((s12 - s13)*TMath::Power(Q2,2)*TMath::Power(s12 - 2*s13,2) + 
         TMath::Power(s12,2)*(-8*s13*TMath::Power(s12,2) + TMath::Power(s12,3) + 11*s12*TMath::Power(s13,2) - 4*TMath::Power(s13,3)) + 
         Q2*s12*(15*s13*TMath::Power(s12,2) - 2*TMath::Power(s12,3) - 25*s12*TMath::Power(s13,2) + 10*TMath::Power(s13,3))) + 
      (s12 - s13)*TMath::Power(s13,2)*TMath::Power(Q2*(s12 - 2*s13) + s12*(-s12 + s13),2)))/fsw2.Re();
  return TMath::Max(0.,ME2);        
}
//____________________________________________________________________________
double Born::Lambda(double a, double b, double c)
{
  return a*a + b*b + c*c - 2*a*b - 2*a*c - 2*b*c;
}
//____________________________________________________________________________
double Born::GetT(double m1, double m2, double m3, double m4, double s, double costh)
{
  //http://edu.itp.phys.ethz.ch/hs10/ppp1/PPP1_2.pdf [Sec 2.2.1]
  double sum = m1*m1+m2*m2+m3*m3+m4*m4;
  return ( (TMath::Sqrt(Lambda(s,m1*m1,m2*m2)*Lambda(s,m3*m3,m4*m4))*costh-(m1*m1-m2*m2)*(m3*m3-m4*m4))/s + sum - s ) /2.;
}