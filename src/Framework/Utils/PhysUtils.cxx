//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 22, 2008 - CA
   That file was added in 2.3.1 - Copied Whitlow R from NuclearUtils and
   added R99118 which was adapted from fortran code sent by V.Tvaskis
 @ Feb 08, 2013 - CA
   Add here the formation zone code so that it can be easily shared between
   the event generation and event reweighting code.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Framework/Utils/PhysUtils.h"
#include "Framework/Messenger/Messenger.h"

//___________________________________________________________________________
double  genie::utils::phys::FormationZone(
   double m, const TLorentzVector & p4, 
   const TVector3 & p3hadr, double ct0 /*in fm*/, double K)
{
// m -> hadon mass (on-shell)
// p -> hadron momentum 4-vector (Lab)
// p3hadr -> hadronic-system momentum 3-vector (Lab)

  TVector3 p3  = p4.Vect();            // hadron's: p (px,py,pz)
  double   m2  = m*m;                  //           m^2
  double   P   = p4.P();               //           |p|
  double   Pt  = p3.Pt(p3hadr);        //           pT
  double   Pt2 = Pt*Pt;                //           pT^2
  double   fz  = P*ct0*m/(m2+K*Pt2);   //           formation zone, in fm

  LOG("PhysUtil", pNOTICE)
      << "Formation zone(|P| = " << P << " GeV, Pt = " << Pt
      << " GeV, ct0 = " << ct0 << " fm, K = " << K << ") = " << fz << " fm";

  return fz;
}
//___________________________________________________________________________
double genie::utils::phys::R99118(double x, double Q2)
{
// Adapted from fortran code sent by V.Tvaskis
// PRL 98, 142301, 2007
//
    double A[3] = { .06723, .46714, 1.89794 };
    double B[3] = { .06347, .57468, -.35342 };
    double C[3] = { .05992, .50885, 2.10807 };
    
    double consq2=2.;  // Joining point of r1990 and re99118
    double R=0;

    if(Q2 >= consq2) {
       double FAC   = 1+12.*(Q2/(1.+Q2))*(.125*.125/(x*x+.125*.125));
       double RLOG  = FAC/TMath::Log(Q2/.04);
       double Q2thr = 5.*TMath::Power(1.-x,5);
       double R_A   = A[0]*RLOG + A[1] / TMath::Sqrt( TMath::Sqrt( TMath::Power(Q2,4)+TMath::Power(A[2],4) ));
       double R_B   = B[0]*RLOG + B[1]/Q2 + B[2]/( TMath::Power(Q2,2) + TMath::Power(.3,2) );
       double R_C   = C[0]*RLOG + C[1]/TMath::Sqrt( TMath::Power(Q2-Q2thr,2) + TMath::Power(C[2],2) );
       R = (R_A+R_B+R_C)/3.;
    }
   
    if(Q2 < consq2) {
       double FAC   = 1+12.*(consq2/(1.+consq2))*(.125*.125/(x*x+.125*.125));
       double RLOG  = FAC/TMath::Log(consq2/.04);
       double Q2thr = 5.*TMath::Power(1.-x,5);
       double R_A   = A[0]*RLOG + A[1]/TMath::Sqrt(TMath::Sqrt( TMath::Power(consq2,4)+TMath::Power(A[2],4) ));
       double R_B   = B[0]*RLOG + B[1]/consq2 + B[2]/( TMath::Power(consq2,2) + TMath::Power(.3,2) );
       double R_C   = C[0]*RLOG + C[1]/TMath::Sqrt( TMath::Power(consq2-Q2thr,2) + TMath::Power(C[2],2) );
       double norm=(R_A+R_B+R_C)/3.;
       R=norm*(1.0-TMath::Exp(-Q2*9.212));
   }
   return R;
}
//___________________________________________________________________________
double genie::utils::phys::RWhitlow(double x, double Q2)
{
// Adapted from NeuGEN's rmodel_mod()
//
// Hugh's comments in original code:
//   from NuTeV code provided by Donna Naples, May 2005
//   added form for R below qsq=.35 GEV**2 from hep-ex/030807
//
// NuTEV comments:
//
//C  Revised to make HTWIST select between more than two different
//C  values of R.
//C
//C     HTWIST = 'F'     ==> RQCD   - WITH LIMIT R > (2MX)**2/Q2
//C     HTWIST = 'T'     ==> RSLAC  - WITH LIMIT R > (2MX)**2/Q2
//C     HTWIST = '0'     ==> R = 0 (NOT TRUE, BUT USEFUL FOR STUDIES)
//C     HTWIST = '2'     ==> R =.2 (NOT TRUE, BUT USEFUL FOR STUDIES)
//C     HTWIST = 'C'     ==> R-CALLAN-GROSS = (2MX)**2/Q2
//C
//C     HTWIST = 'W'     ==> R PARAMETERIZATION FROM WHITLOW'S THESIS
//C     HTWIST = '+'     ==> R PARAMETRIZATION FROM WHITLOW +15%
//C     HTWIST = 'P'     ==> R PARAMETRIZATION FROM WHITLOW +0.03
//C     HTWIST = 'M'     ==> R PARAMETRIZATION FROM WHITLOW -0.03
//C
//C  22-DEC-92 WGS
//
  const double C2 = TMath::Power(0.125, 2);
  const double B1 =  0.0635;
  const double B2 =  0.5747;
  const double B3 = -0.3534;

  double Q2R   = TMath::Max(Q2, 0.35);
  double ss    = TMath::Log(Q2R/.04);

  double x2    = TMath::Power(x,   2.);
  double Q4R   = TMath::Power(Q2R, 2.);
  double Q4    = TMath::Power(Q2,  2.);

  double theta = 1. + (12.*Q2R/(Q2R+1.)) * (C2/(C2+x2));
  double R     = (B1/ss)*theta + B2/Q2R + B3/(Q4R+.09);

  R = TMath::Max(R,0.);

  if(Q2 < 0.35)  {
      R *= ( 3.207*(Q2/(Q4+1.)) );
  }
  return R;
}
//___________________________________________________________________________
/*
void genie::utils::phys::ExtractStructFunc (
        double x, double Q2, double d2sig_dxdy[3], 
        double & F1, double & F2, double & xF3)
{
// Solve the system of equations:
//   (Sigma) = (A) x (SF) => SF = (A)^-1 x (Sigma)
//

  const double sign = 1;
  const double M    = kNucleonMass;

  TMatrixD A     (3,3); // (row,col) 
  TMatrixD Sigma (3,1);

  A(0,0) = x * TMath::Power(y[0],2.);
  A(1,0) = x * TMath::Power(y[1],2.);
  A(2,0) = x * TMath::Power(y[2],2.);

  A(0,1) = 1. - y[0] - 0.5*M*x*y[0]/E[0];
  A(1,1) = 1. - y[1] - 0.5*M*x*y[1]/E[1];
  A(2,1) = 1. - y[2] - 0.5*M*x*y[2]/E[2];

  A(0,2) = sign * y[0] * (1.-0.5*y[0]);
  A(1,2) = sign * y[1] * (1.-0.5*y[1]);
  A(2,2) = sign * y[2] * (1.-0.5*y[2]);

  Sigma(0,0) = d2sig_dxdy[0] / (kGF2*M*E[0]/kPi);
  Sigma(1,0) = d2sig_dxdy[1] / (kGF2*M*E[1]/kPi);
  Sigma(2,0) = d2sig_dxdy[2] / (kGF2*M*E[2]/kPi);

  A.Invert();

  TMatrixD SF = A*Sigma;

  F1  = SF(0,0);
  F2  = SF(1,0);
  xF3 = SF(2,0);
}
//___________________________________________________________________________
*/
