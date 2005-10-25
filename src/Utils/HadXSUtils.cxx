//____________________________________________________________________________
/*!

\namespace  genie::utils::hadxs

\brief      Simple functions and data for computing hadron interaction xsecs

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    March 11, 2004

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Utils/HadXSUtils.h"

#include <TMath.h>

using namespace genie::constants;

namespace genie {
namespace utils {
namespace hadxs {

 // The inelastic pion-nucleon cross section used in Rein-Seghal coherent pi0
 // production xsec: D.Rein and L.M.Seghal,Nucl.Phys.B223:29-144 (1983).
 // The data used here are coming from CERN-HERA 79-01, 1979, 'Compilation of
 // cross sections I - pi- and pi+ induced reactions'. Also, look at the
 // Fig.3 in Rein-Seghal's paper.
 // However, the actual values I am using are copied from Hugh Gallagher's
 // NeuGEN inel() function which is adapted here.
 const int    kInelNDataPoints = 60;
 const double kInelMinLog10P   = -0.975;
 const double kIneldLog10P     =  0.05;
 const double kInelSig[kInelNDataPoints] = {
         1.596,   3.192,   5.692,   5.596,   3.788,   6.528,
        22.931,  43.462,  55.580,  36.761,  19.754,  12.588,
        12.914,  12.500,  13.707,  15.082,  19.671,  16.860,
        21.708,  29.128,  21.752,  22.444,  23.698,  23.847,
        23.067,  25.336,  25.366,  25.273,  24.646,  24.003,
        23.636,  23.615,  23.029,  22.667,  22.434,  21.901,
        21.763,  22.235,  20.177,  21.707,  20.827,  21.102,
        21.028,  21.155,  20.932,  20.577,  20.865,  21.122,
        21.193,  21.081,  20.611,  20.788,  20.591,  20.514,
        20.796,  20.813,  20.425,  20.460,  20.495,  20.530
 };

 // Total pion nucleon cross section data
 // Data from Hugh Gallagher's NeuGEN total() function.
 const int    kTotNDataPoints = 60;
 const double kTotMinLog10P   = -0.975;
 const double kTotdLog10P     =  0.05;
 const double kTotSig[kTotNDataPoints] = {
         3.252,   6.504,  12.316,  18.314,  22.600,  31.435,
        53.933,  84.872, 102.626,  87.084,  60.234,  39.922,
        32.804,  28.935,  27.952,  29.439,  36.824,  31.814,
        35.152,  50.062,  39.079,  35.741,  37.390,  35.575,
        34.043,  34.363,  34.171,  32.990,  32.110,  31.316,
        30.621,  29.918,  29.134,  28.461,  27.985,  27.208,
        26.749,  27.111,  25.047,  26.357,  25.393,  25.777,
        25.467,  25.305,  25.008,  24.814,  24.662,  24.481,
        24.453,  24.336,  24.099,  24.098,  24.010,  23.908,
        23.992,  24.058,  23.805,  23.808,  23.811,  23.815
 };

}
}
}
//____________________________________________________________________________
double genie::utils::hadxs::InelasticPionNucleonXSec(double Epion)
{
// Returns the inelastic pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN inel() function.

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass_2) );

  double log10P = 0;
  if(P>0) log10P = TMath::Log10(P);
  else return 0.;

  int N = (int) ((log10P - kInelMinLog10P)/kIneldLog10P);

  if      (N<0)                  return (P/0.1059)*kInelSig[0];
  else if (N>kInelNDataPoints-2) return kInelSig[kInelNDataPoints-1];
  else {
   double dsig = kInelSig[N+1]-kInelSig[N];
   return (log10P-kInelMinLog10P-N*kIneldLog10P)*dsig/kIneldLog10P;
  }
}
//____________________________________________________________________________
double genie::utils::hadxs::TotalPionNucleonXSec(double Epion)
{
// Returns the total pion-nucleon cross section.
// C++ adaptation of Hugh Gallagher's NeuGEN total() function.

  double Epion2 = TMath::Power(Epion,2);
  double P      = TMath::Sqrt( TMath::Max(0.,Epion2-kPionMass_2) );

  double log10P = 0;
  if(P>0) log10P = TMath::Log10(P);
  else    return 0.;

  int N = (int) ((log10P - kTotMinLog10P)/kTotdLog10P);

  if      (N<0)                  return (P/0.1059)*kTotSig[0];
  else if (N>kInelNDataPoints-2) return kTotSig[kInelNDataPoints-1];
  else {
   double dsig = kTotSig[N+1]-kTotSig[N];
   return (log10P-kTotMinLog10P-N*kTotdLog10P)*dsig/kTotdLog10P;
  }
}
//____________________________________________________________________________

