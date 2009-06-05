//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         May 30, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 30, 2009 - CA
   This file was added in v2.5.1.

*/
//____________________________________________________________________________

#include <TMatrixD.h>

#include "ValidationTools/StructFunc/ExtractStructFunc.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::vld_structfunc;

//____________________________________________________________________________
void genie::vld_structfunc::ExtractStructFunc(
         double x, double Q2,  int lepton, int nucleon, 
         double & F1, double & F2, double & xF3)
{
// Extract F1, F2, xF3 from a "black box" cross section model as described in 
// H.Gallagher, Nucl.Phys.Proc.Suppl.159:229-234,2006
//
  const double sign = 1;
  const double M    = kNucleonMass;

  const int N = 3;

  double E[N]          = {4., 10., 40.};
  double y[N]          = {0.,  0.,  0.};
  double d2sig_dxdy[N] = {0.,  0.,  0.};

  // E,y pairs for constant x,Q2 ( Q^2/2Mx = y*E)
  for(int i=0; i<N; i++) {
     y[i] = Q2/(2*M*x) / E[i];
  }

  // Calculate the total d2sigma/dxdy cross section for each {x,Q2,y[i],E[i]}
  for(int i=0; i<N; i++) {

     // add code here to calculate d2sigma/dxdy
     //
     d2sig_dxdy[i] = 1;
  }
  
  //
  // Solve the system of equations:
  //   (Sigma) = (A) x (SF) => SF = (A)^-1 x (Sigma)
  //
  TMatrixD A     (3,3); /* (row,col) */
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
//____________________________________________________________________________
