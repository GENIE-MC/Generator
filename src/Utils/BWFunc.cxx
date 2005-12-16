//____________________________________________________________________________
/*!

\file     breit_wigner_func

\brief    Breit Wigner functions

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Utils/BWFunc.h"
#include "Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

//______________________________________________________________________
double genie::utils::bwfunc::BreitWignerL(
               double W, int L, double mass, double width0, double norm)
{
//Inputs:
// - W:      Invariant mass (GeV)
// - L:      Resonance orbital angular momentum
// - mass:   Resonance mass (GeV)
// - width0: Resonance width
// - norm:   Breit Wigner norm

  //-- sanity checks
  assert(mass   >  0);
  assert(width0 >  0);
  assert(norm   >  0);
  assert(W      >  0);
  assert(L      >= 0);

  //-- auxiliary parameters
  double mN    = kNucleonMass;
  double mPi   = kPionMass;
  double m_2   = TMath::Power(mass, 2);
  double mN_2  = TMath::Power(mN,   2);
  double mPi_2 = TMath::Power(mPi,  2);
  double W_2   = TMath::Power(W,    2);

  //-- calculate the L-dependent resonance width
  double qpW_2 = ( TMath::Power(W_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  double qpM_2 = ( TMath::Power(m_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  if(qpW_2 < 0) qpW_2 = 0;
  if(qpM_2 < 0) qpM_2 = 0;
  double qpW   = TMath::Sqrt(qpW_2) / (2*W);
  double qpM   = TMath::Sqrt(qpM_2) / (2*mass);
  double width = width0 * TMath::Power( qpW/qpM, 2*L+1 );

  //-- calculate the Breit Wigner function for the input W
  double width_2 = TMath::Power( width,  2);
  double W_m_2   = TMath::Power( W-mass, 2);

  double bw = (0.5/kPi) * (width/norm) / (W_m_2 + 0.25*width_2);
  return bw;
}
//______________________________________________________________________
double genie::utils::bwfunc::BreitWigner(
                       double W, double mass, double width, double norm)
{
//Inputs:
// - W:     Invariant mass (GeV)
// - mass:  Resonance mass (GeV)
// - width: Resonance width
// - norm:  Breit Wigner norm

  //-- sanity checks
  assert(mass  >  0);
  assert(width >  0);
  assert(norm  >  0);
  assert(W     >  0);

  //-- auxiliary parameters
  double width_2 = TMath::Power( width,  2);
  double W_m_2   = TMath::Power( W-mass, 2);

  //-- calculate the Breit Wigner function for the input W
  double bw = (0.5/kPi) * (width/norm) / (W_m_2 + 0.25*width_2);
  return bw;
}
//______________________________________________________________________

