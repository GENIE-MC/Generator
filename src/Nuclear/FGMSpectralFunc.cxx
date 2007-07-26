//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - June 20, 2004

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
bool FGMSpectralFunc::GenerateNucleon(const Target & /*target*/) const
{
  return false;
}
//____________________________________________________________________________
double FGMSpectralFunc::Prob(
          double /*p*/, double /*w*/, const Target & /*target*/) const
{
  return 0;
}
//____________________________________________________________________________
