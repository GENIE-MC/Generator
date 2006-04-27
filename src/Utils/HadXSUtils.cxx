//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 11, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Utils/HadXSUtils.h"

#include <TMath.h>

using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::hadxs::InelasticPionNucleonXSec(double Epion)
{
// Returns the inelastic pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN inel() function.
// Actually, the following is a simple data interpolation:

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass2) );

  double log10P = 0;
  if(P>0) log10P = TMath::Log10(P);
  else return 0.;

  int N = (int) ((log10P - kInelMinLog10P)/kIneldLog10P);

  double log10Pn = kInelMinLog10P +  N * kIneldLog10P;

  if      (N<0)                  return (P/0.1059)*kInelSig[0];
  else if (N>kInelNDataPoints-2) return kInelSig[kInelNDataPoints-1];
  else {
   double d  = (kInelSig[N+1]-kInelSig[N])/kIneldLog10P;
   double xs = kInelSig[N] + d * (log10P-log10Pn);
   return xs;
  }
}
//____________________________________________________________________________
double genie::utils::hadxs::TotalPionNucleonXSec(double Epion)
{
// Returns the total pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN total() function.
// Actually, the following is a simple data interpolation:

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass2) );

  double log10P = 0;
  if(P>0) log10P = TMath::Log10(P);
  else    return 0.;

  int N = (int) ((log10P - kTotMinLog10P)/kTotdLog10P);
  double log10Pn = kTotMinLog10P +  N * kTotdLog10P;

  if      (N<0)                  return (P/0.1059)*kTotSig[0];
  else if (N>kInelNDataPoints-2) return kTotSig[kInelNDataPoints-1];
  else {
   double d  = (kTotSig[N+1]-kTotSig[N])/kTotdLog10P;
   double xs = kTotSig[N] + d * (log10P-log10Pn);
   return xs;
  }
}
//____________________________________________________________________________

