//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - June 20, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Nuclear/FGMSpectralFunc.h"

using namespace genie;

//____________________________________________________________________________
FGMSpectralFunc::FGMSpectralFunc() :
NuclearModelI("genie::FGMSpectralFunc")
{

}
//____________________________________________________________________________
FGMSpectralFunc::FGMSpectralFunc(string config) :
NuclearModelI("genie::FGMSpectralFunc", config)
{

}
//____________________________________________________________________________
FGMSpectralFunc::~FGMSpectralFunc()
{

}
//____________________________________________________________________________
bool FGMSpectralFunc::GenerateNucleon(const Target & t) const
{
  return false;
}
//____________________________________________________________________________
double FGMSpectralFunc::Prob(double p, double w, const Target & t) const
{
  return 0;
}
//____________________________________________________________________________
