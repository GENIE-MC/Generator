//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - June 15, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cmath>

#include "Fragmentation/FragmentationFunctions.h"

//___________________________________________________________________________
double genie::utils::frgmfunc::collins_spiller_func(double* x, double* par)
{
// par[0] = N
// par[1] = epsilon

  double z = x[0];

  double D = par[0] * ( (1.-z)/z + par[1]*(2.-z)/(1.-z) ) *
                         pow(1+z, 2.) * pow(1. - 1./z - par[1]/(1.-z), -2.);
  return D;
}
//___________________________________________________________________________
double genie::utils::frgmfunc::peterson_func(double* x, double* par)
{
// par[0] = N
// par[1] = epsilon

  double z = x[0];

  double D = par[0] / ( z * pow(1. - 1./z - par[1]/(1.-z), 2) );

  return D;
}
//___________________________________________________________________________

