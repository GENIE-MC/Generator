//____________________________________________________________________________
/*!

\file     fragmentation_functions.cxx

\brief    Fragmentation functions

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004
 
*/
//____________________________________________________________________________

#include <cmath>

#include "fragmentation_functions.h"

using namespace genie;

//___________________________________________________________________________
double genie::collins_spiller_fragmentation_function(double * x, double * par)
{
// par[0] = N
// par[1] = epsilon

  double z = x[0];

  double D = par[0] * ( (1.-z)/z + par[1]*(2.-z)/(1.-z) ) *
                         pow(1+z, 2.) * pow(1. - 1./z - par[1]/(1.-z), -2.);
  return D;
}
//___________________________________________________________________________
double genie::peterson_fragmentation_function(double * x, double * par)
{
// par[0] = N
// par[1] = epsilon

  double z = x[0];

  double D = par[0] / ( z * pow(1. - 1./z - par[1]/(1.-z), 2) );

  return D;
}
//___________________________________________________________________________

