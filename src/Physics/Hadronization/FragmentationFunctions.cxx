//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cmath>

#include "Physics/Hadronization/FragmentationFunctions.h"

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
