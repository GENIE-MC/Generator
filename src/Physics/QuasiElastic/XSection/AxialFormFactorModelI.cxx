//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Aaron Meyer <asmeyer2012 \at uchicago.edu>

 based off ELFormFactorsModelI by
 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

using namespace genie;

//____________________________________________________________________________
AxialFormFactorModelI::AxialFormFactorModelI() :
Algorithm()
{

}
//____________________________________________________________________________
AxialFormFactorModelI::AxialFormFactorModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
AxialFormFactorModelI::AxialFormFactorModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
AxialFormFactorModelI::~AxialFormFactorModelI()
{

}
//____________________________________________________________________________
