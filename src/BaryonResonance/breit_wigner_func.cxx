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

#include "BaryonResonance/breit_wigner_func.h"
#include "Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

//______________________________________________________________________
double genie::breit_wigner_L(double W, const Registry & res_config)
{
  //-- get mass, width, norm and orbital angular momentum from config.
  
  assert(
      res_config.Exists("Res-Orb-Angular-Mom") &&
      res_config.Exists("Res-Mass")            &&
      res_config.Exists("Res-Width")           &&
      res_config.Exists("Breit-Wigner-Norm")
  );

  int    L       =  res_config.GetInt    ( "Res-Orb-Angular-Mom");
  double mass    =  res_config.GetDouble ( "Res-Mass"           );
  double width0  =  res_config.GetDouble ( "Res-Width"          );
  double bw_norm =  res_config.GetDouble ( "Breit-Wigner-Norm"  );

  //-- sanity checks

  assert(mass    >  0);
  assert(width0  >  0);
  assert(bw_norm >  0);
  assert(W       >  0);
  assert(L       >= 0);
    
  //-- auxiliary parameters
    
  double mN  = kNucleonMass;
  double mPi = kPionMass;

  double m_2   = pow(mass, 2);
  double mN_2  = pow(mN,   2);
  double mPi_2 = pow(mPi,  2);
  double W_2   = pow(W,    2);
  
  //-- calculate the L-dependent resonance width

  double qpW_2 = ( pow(W_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  double qpM_2 = ( pow(m_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  
  if(qpW_2 < 0) qpW_2 = 0;
  if(qpM_2 < 0) qpM_2 = 0;
  
  double qpW   = sqrt(qpW_2) / (2*W);
  double qpM   = sqrt(qpM_2) / (2*mass);

  double width  = width0 * pow( qpW/qpM, 2*L+1 );
  
  //-- calculate the Breit Wigner function for the input W

  double width_2 = pow( width,  2);
  double W_m_2   = pow( W-mass, 2);
  
  double bw = (0.5/kPi) * (width/bw_norm) / (W_m_2 + 0.25*width_2);

  return bw;
}
//______________________________________________________________________
double genie::breit_wigner(double W, const Registry & res_config)
{
  //-- get mass, width, norm and orbital angular momentum from config.

  assert(
      res_config.Exists("Res-Mass")          &&
      res_config.Exists("Res-Width")         &&
      res_config.Exists("Breit-Wigner-Norm")
  );

  double mass    =  res_config.GetDouble ( "Res-Mass"           );
  double width   =  res_config.GetDouble ( "Res-Width"          );
  double bw_norm =  res_config.GetDouble ( "Breit-Wigner-Norm"  );

  //-- sanity checks

  assert(mass    >  0);
  assert(width   >  0);
  assert(bw_norm >  0);
  assert(W       >  0);

  //-- auxiliary parameters

  double width_2 = pow( width,  2);
  double W_m_2   = pow( W-mass, 2);

  //-- calculate the Breit Wigner function for the input W

  double bw = (0.5/kPi) * (width/bw_norm) / (W_m_2 + 0.25*width_2);

  return bw;
}
//______________________________________________________________________

