//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
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
  
  fgae = -1./2. + 2.*fsw2;
  fgbe = -1./2.;
  fgav = 1./2.;

}
//____________________________________________________________________________
Born::~Born()
{

/*
Make sure the p3 is always the charged lepton:

nu (p1) + lp (p2) -> lp (p3) + nu (p4)

  s-channel       t-channel        u-channel
1 \      / 3    1 --------- 3    1 ------\ / 3
   ------             |                |  X
2 /      \ 4    2 --------- 4    2 ------/ \ 4
*/

}
//____________________________________________________________________________
double Born::PXSecCCR(double s, double t, double mlin, double mlout)
{

  TComplex prop = falpha/fsw2/(s-fmw2c);

  return (t-mlout*mlout)*(t-mlin*mlin) * prop.Rho2();

}
//____________________________________________________________________________
double Born::PXSecCCV(double s, double t, double mlin, double mlout)
{

  TComplex prop = falpha/fsw2/(t-fmw2c);

  return (s-mlout*mlout)*(s-mlin*mlin) * prop.Rho2();

}
//____________________________________________________________________________
double Born::PXSecCCRNC(double s, double t, double mlin, double mlout)
{

  double u = GetU(mlin,mlout,s,t);
  
  TComplex a = fgav*(fgae-fgbe)/(u-fmz2c)/fcw2/fsw2;
  TComplex b = fgav*(fgae+fgbe)/(u-fmz2c)/fcw2/fsw2 + 1./(s-fmw2c)/fsw2;
  return falpha.Rho2() * ( (t-mlout*mlout)*(t-mlin*mlin)*b.Rho2() + (s-mlout*mlout)*(s-mlin*mlin)*a.Rho2() );

}
//____________________________________________________________________________
double Born::PXSecCCVNC(double s, double t, double mlin, double mlout)
{

  double u = GetU(mlin,mlout,s,t);
  
  TComplex a = fgav*(fgae+fgbe)/(u-fmz2c)/fcw2/fsw2 + 1./(t-fmw2c)/fsw2;
  TComplex b = fgav*(fgae-fgbe)/(u-fmz2c)/fcw2/fsw2;
  return falpha.Rho2() * ( (t-mlout*mlout)*(t-mlin*mlin)*b.Rho2() + (s-mlout*mlout)*(s-mlin*mlin)*a.Rho2() );

}
//____________________________________________________________________________
double Born::PXSecNCVnu(double s, double t, double mlin, double mlout)
{

  double u = GetU(mlin,mlout,s,t);
  
  TComplex a = fgav*(fgae+fgbe)/(u-fmz2c)/fcw2/fsw2;
  TComplex b = fgav*(fgae-fgbe)/(u-fmz2c)/fcw2/fsw2;
  return falpha.Rho2() * ( (t-mlout*mlout)*(t-mlin*mlin)*b.Rho2() + (s-mlout*mlout)*(s-mlin*mlin)*a.Rho2() );

}
//____________________________________________________________________________
double Born::PXSecNCVnubar(double s, double t, double mlin, double mlout)
{

  double u = GetU(mlin,mlout,s,t);
  
  TComplex a = fgav*(fgae-fgbe)/(u-fmz2c)/fcw2/fsw2;
  TComplex b = fgav*(fgae+fgbe)/(u-fmz2c)/fcw2;
  return falpha.Rho2() * ( (t-mlout*mlout)*(t-mlin*mlin)*b.Rho2()/fsw2.Rho2() + (s-mlout*mlout)*(s-mlin*mlin)*a.Rho2() );

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
double Born::PXSecPhoton_T(double s12, double s13, double Q2, double ml2)
{
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
double Born::PXSecPhoton_L(double s12, double s13, double Q2, double ml2)
{
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
double Born::GetS(double mlin, double Enuin)
{
  return 2. * mlin * Enuin + mlin*mlin;
}
//____________________________________________________________________________
double Born::GetT(double mlin, double mlout, double s, double costhCM)
{
  //http://edu.itp.phys.ethz.ch/hs10/ppp1/PPP1_2.pdf [Sec 2.2.1]
  double sum = mlin*mlin+mlout*mlout;
  return ( (TMath::Sqrt(Lambda(s,0.,mlin*mlin)*Lambda(s,mlout*mlout,0.))*costhCM+mlin*mlin*mlout*mlout)/s + sum - s ) /2.;
}
//____________________________________________________________________________
double Born::GetU(double mlin, double mlout, double s, double t)
{
    return mlin*mlin+mlout*mlout-s-t;
}
//____________________________________________________________________________
bool Born::IsInPhaseSpace(double mlin, double mlout, double Enuin, double Enuout)
{

  //https://arxiv.org/pdf/2007.14426.pdf [section 2.2]
  double frac = Enuout/Enuin;
  if      ( frac < mlin/(mlin+2.*Enuin)+(mlout*mlout-mlin*mlin)/2./Enuin/(mlin+2.*Enuin)  ) return false;
  else if ( frac > 1.-(mlout*mlout-mlin*mlin)/2./Enuin/mlin ) return false;

  return true;

}




