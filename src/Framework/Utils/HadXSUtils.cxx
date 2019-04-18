//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - March 11, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jul 18, 2013 - Daniel Scully
 Fixed indexing bug in InelasticPionNucleonXSec and TotalPionNucleonXSec
 @ Jul 18, 2013 - Daniel Scully
 Added pion-nucleon cross-sections from C. Berger, provided via D. Cherdack
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Utils/HadXSUtils.h"

using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::hadxs::InelasticPionNucleonXSec(double Epion,
    bool isChargedPion)
{
// Returns the interpolated inelastic pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN inel() function

  double mpi  = kPionMass;
  double mpi2 = kPionMass2;
  if (!isChargedPion) {
    mpi  = kPi0Mass;
    mpi2 = mpi * mpi;
  }
  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-mpi2) );

  if(P<=0) return 0;

  double log10P  = TMath::Log10(P);
  int    N = (int) ((log10P - kInelMinLog10P)/kIneldLog10P) + 1;

  double xs=0.;
  if ((log10P - kInelMinLog10P) < 0.0) xs = (P/0.1059)*kInelSig[0];
  else if (N>kInelNDataPoints-2) xs = kInelSig[kInelNDataPoints-1];
  else {
   double log10Pn = kInelMinLog10P +  (N-1) * kIneldLog10P;
   double delta   = (kInelSig[N]-kInelSig[N-1])/kIneldLog10P;
   xs = kInelSig[N-1] + delta * (log10P-log10Pn);
  }
  return (xs * units::mb);
}
//____________________________________________________________________________
double genie::utils::hadxs::TotalPionNucleonXSec(double Epion,
    bool isChargedPion)
{
// Returns the interpolated total pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN total() function.

  double mpi  = kPionMass;
  double mpi2 = kPionMass2;
  if (!isChargedPion) {
    mpi  = kPi0Mass;
    mpi2 = mpi * mpi;
  }
  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-mpi2) );

  if(P<=0) return 0;

  double log10P  = TMath::Log10(P);
  int    N = (int) ((log10P - kInelMinLog10P)/kIneldLog10P) + 1;

  double xs=0.;
  if ((log10P - kInelMinLog10P) < 0.0) xs = (P/0.1059)*kTotSig[0];
  else if (N>kInelNDataPoints-2) xs = kTotSig[kInelNDataPoints-1];
  else {
   double log10Pn = kTotMinLog10P +  (N-1) * kTotdLog10P;
   double delta   = (kTotSig[N]-kTotSig[N-1])/kTotdLog10P;
   xs = kTotSig[N-1] + delta * (log10P-log10Pn);
  }
  return (xs * units::mb);
}
//____________________________________________________________________________
double genie::utils::hadxs::berger::InelasticPionNucleonXSec(double Epion, 
    bool isChargedPion)
{
  const double total = PionNucleonXSec(Epion, true, isChargedPion);
  const double elastic = PionNucleonXSec(Epion, false, isChargedPion);
  return (total - elastic);
}
//____________________________________________________________________________
double genie::utils::hadxs::berger::TotalPionNucleonXSec(double Epion, 
    bool isChargedPion)
{
  return PionNucleonXSec(Epion, true, isChargedPion);
}
//____________________________________________________________________________
double genie::utils::hadxs::berger::PionNucleonXSec(double Epion, bool get_total, 
    bool isChargedPion)
{
  // Convert inputs from Genie to those expected by Berger's code:

  double mpi  = kPionMass;
  double mpi2 = kPionMass2;
  if (!isChargedPion) {
    mpi  = kPi0Mass;
    mpi2 = mpi * mpi;
  }
  
  double Epion2 = TMath::Power(Epion,2);
  double ppi = TMath::Sqrt( TMath::Max(0., Epion2 - mpi2) );
  
  if( ppi <= 0.0 ) return 0.0;
  
  int out;
  if( get_total ) out = 0;  // Total pion-nucleon cross-section
  else out = 1;             // Elastic pion-nucleon cross-section
  
  const double M_pi = mpi;        // TODO: used to be kPionMass prior to charge checks
  const double M_p = kProtonMass; // TODO: should be kNucleonMass  ??
  const double pi = kPi;
  
  // Now this is the Berger's code...
  
  double afit=1.0;
  double afit2=1.9;
  double afit3=0.27;
  double afit4=0.34;
  double afit5=0.75;
  double afit6=1.7;

  double epi = TMath::Sqrt(M_pi*M_pi + ppi*ppi);
  double sx  = M_p*M_p + M_pi*M_pi + 2.0 * epi * M_p;
  double Wx  = TMath::Sqrt(sx);

  double s12x   = TMath::Sqrt((sx - TMath::Power((M_pi+M_p),2)) * (sx - (M_pi-M_p)*(M_pi-M_p)));
  double ppistx = s12x/2.0/Wx;

  double Wdel  = 1.232;
  double arg1  = (Wdel*Wdel - (M_p + M_pi)*(M_p + M_pi))*(Wdel*Wdel - (M_p -M_pi)*(M_p -M_pi));
  if(arg1<=0) arg1 = 0.0;
  double gamx  = 0.15 * TMath::Power((6.0*ppistx),3) / (1.0 + (6*ppistx)*(6*ppistx));
  double bwnrx = 1.0 / (4.0*(Wx - Wdel)*(Wx - Wdel) + gamx*gamx);
  double f1x   = afit * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gamx*gamx * bwnrx;
  double f4x   = afit4 * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gamx*gamx * bwnrx;

  double Wres2  = 1.98;
  double Gres2  = 0.6;
  double arg2   = (Wres2*Wres2 - (M_p+M_pi)*(M_p+M_pi)) * (Wres2*Wres2 - (M_p - M_pi)*(M_p - M_pi));
  if(arg2<=0) arg2 = 0.0;
  double pp2    = TMath::Sqrt(arg2)/2.0/Wres2;
  double gam2x  = Gres2 * TMath::Power((ppistx/pp2),3);
  double bwnr2x = 1.0 / (4.0 * (Wx-Wres2)*(Wx-Wres2) + gam2x*gam2x);
  double f2x    = afit2 * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gam2x*gam2x * bwnr2x;

  double Wres3  = 1.7;
  double Gres3  = 0.27;
  double arg3   = (Wres3*Wres3 - (M_p+M_pi)*(M_p+M_pi)) * (Wres3*Wres3 - (M_p - M_pi)*(M_p - M_pi));
  if(arg3<=0) arg3 = 0.0;
  double pp3    = TMath::Sqrt(arg3)/2.0/Wres3;
  double gam3x  = Gres3 * TMath::Power((ppistx/pp3),3);
  double bwnr3x = 1.0 / (4.0 * (Wx-Wres3)*(Wx-Wres3) + gam3x*gam3x);
  double f3x    = afit3 * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gam3x*gam3x * bwnr3x;

  double Wres5  = 1.52;
  double Gres5  = 0.10;
  double arg5   = (Wres5*Wres5 - (M_p+M_pi)*(M_p+M_pi)) * (Wres5*Wres5 - (M_p - M_pi)*(M_p - M_pi));
  if(arg5<=0) arg5 = 0.0;
  double pp5    = TMath::Sqrt(arg5)/2.0/Wres5;
  double gam5x  = Gres5 * (ppistx/pp5);
  double bwnr5x = 1.0 / (4.0 * (Wx-Wres5)*(Wx-Wres5) + gam5x*gam5x);
  double f5x    = afit5 * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gam5x*gam5x * bwnr5x;

  double Wres6  = 1.69;
  double Gres6  = 0.140;
  double arg6   = (Wres6*Wres6 - (M_p+M_pi)*(M_p+M_pi)) * (Wres6*Wres6 - (M_p - M_pi)*(M_p - M_pi));
  if(arg6<=0) arg6 = 0.0;
  double pp6    = TMath::Sqrt(arg6)/2.0/Wres6;
  double gam6x  = Gres6 * (ppistx/pp6);
  double bwnr6x = 1.0 / (4.0 * (Wx-Wres6)*(Wx-Wres6) + gam6x*gam6x);
  double f6x    = afit6 * 0.3892 * 8.0 * pi / (ppistx*ppistx) * gam6x*gam6x * bwnr6x;

  double output = -9999.99;
  if(out==0){
    double bgplusx  = 0.834 * ppi + 1.77 * ppi*ppi;
    double sigplusx = f1x + 0.86 * f2x + 0.7 * f3x + bgplusx;
    if(ppi>1.9) sigplusx = 20.1 + 15.4 / TMath::Sqrt(ppi);
    double bgminx  = 17.8 * ppi + 6.22 * ppi*ppi - 3.1 * ppi*ppi*ppi;
    double sigminx = 0.94 * f4x + 0.65 *f5x + 0.61 * f6x + bgminx;
    if(ppi>1.2) sigminx = 21.64 + 17.29 / TMath::Sqrt(ppi);
    double sigma_t = (sigplusx + sigminx) / 2.0;
    output = sigma_t;
  }else if(out==1){
    double sgpelx = f1x + 0.411 * f2x;
    if(ppi>1.9) sgpelx = 2.38 + 17.55 / ppi;
    double sgmelx = 0.275 * f4x + 0.268 * f5x + 0.281 * f6x + 11.83 * ppi - 3.39 * ppi*ppi;
    if(ppi>1.5) sgmelx = 3.4 + 11.35 / ppi;
    double sigma_el = (sgpelx + sgmelx) / 2.0;
    output = sigma_el;
  }
  
  // Berger's code over, convert to Genie units and return
  return (output * units::mb);
}
//____________________________________________________________________________
// Berger code for calculating the pion-Carbon xsec. If A is not 12 the restults are extrapolated to a nucleus of size A.
int genie::utils::hadxs::berger::PionNucleusXSec(double tpi, double ppistar, double t_new, double A, double &tpilow, double &siglow, double &tpihigh, double &sighigh){

  //Berger code for Tpi<=1.0
  //Returns the entire dsigma(pi + N -> pi + N)/dt term based on pi-Carbon scattering data
  double binedges[13] = {0.000, 0.076, 0.080, 0.100, 0.148, 0.162, 0.226, 0.486, 0.584, 0.662, 0.766, 0.870, 1.000};
  double parones[13]  = {0.0, 11600.0, 14700.0, 18300.0, 21300.0, 22400.0, 16400.0, 5730.0, 4610.0, 4570.0, 4930.0, 5140.0, 5140.0};
  double partwos[13]  = {0.0, 116.0, 109.0, 89.8, 91.0, 89.2, 80.8, 54.6, 55.2, 58.4, 60.5, 62.2, 62.2};
  double factors[13]  = {0.0, 0.1612, 0.1662, 0.1906, 0.2452, 0.2604, 0.3273, 0.5784, 0.6682, 0.7384, 0.8304, 0.9206, 0.9206};

  if(tpi>binedges[12]) return 1;
  int btu = 1;
  while(true){
    if(tpi>binedges[btu]) btu++;
    else                  break;
  }
  btu--;
  if(btu<0) btu = 0;

  tpilow  = binedges[btu];
  tpihigh = binedges[btu + 1];

  double dsigdzlow = -9999.99;
  if(btu==0) dsigdzlow = 0.0;
  else       dsigdzlow = 2.0 * factors[btu]*factors[btu]
    * parones[btu]*TMath::Power(A/12.0,1.3333333)
      * TMath::Exp( -1.0 * partwos[btu]*TMath::Power(A/12.0,0.6666666) * t_new * factors[btu]*factors[btu] / (ppistar*ppistar) );
  siglow = dsigdzlow;

  double dsigdzhigh = -9999.99;
  btu++;
  if(btu>12) btu        = 12;
  if(btu==0) dsigdzhigh = 0.0;
  else       dsigdzhigh = 2.0 * factors[btu]*factors[btu]
    * parones[btu]*TMath::Power(A/12.0,1.3333333)
      * TMath::Exp( -1.0 * partwos[btu]*TMath::Power(A/12.0,0.6666666) * t_new * factors[btu]*factors[btu] / (ppistar*ppistar) );
  sighigh = dsigdzhigh;

  return 0;
}
//_____________________________________________________________________________
