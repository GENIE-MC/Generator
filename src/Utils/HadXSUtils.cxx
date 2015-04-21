//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - March 11, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jul 18, 2013 - Daniel Scully
 Fixed indexing bug in InelasticPionNucleonXSec and TotalPionNucleonXSec

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Utils/HadXSUtils.h"

using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::hadxs::InelasticPionNucleonXSec(double Epion)
{
// Returns the interpolated inelastic pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN inel() function

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass2) );

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
double genie::utils::hadxs::TotalPionNucleonXSec(double Epion)
{
// Returns the interpolated total pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN total() function.

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass2) );

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

