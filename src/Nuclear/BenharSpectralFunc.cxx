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

#include "Nuclear/BenharSpectralFunc.h"

using namespace genie;

//____________________________________________________________________________
BenharSpectralFunc::BenharSpectralFunc() :
NuclearModelI("genie::BenharSpectralFunc")
{

}
//____________________________________________________________________________
BenharSpectralFunc::BenharSpectralFunc(string config) :
NuclearModelI("genie::BenharSpectralFunc", config)
{

}
//____________________________________________________________________________
BenharSpectralFunc::~BenharSpectralFunc()
{

}
//____________________________________________________________________________
bool BenharSpectralFunc::GenerateNucleon(const Target & /*target*/) const
{
  return false;
}
//____________________________________________________________________________
double BenharSpectralFunc::Prob(
             double /*p*/, double /*w*/, const Target & /*target*/) const
{
  return 0;
}
//____________________________________________________________________________
