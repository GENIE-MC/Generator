//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

using namespace genie;

//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI() :
Algorithm()
{

}
//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
ELFormFactorsModelI::~ELFormFactorsModelI()
{

}
//____________________________________________________________________________
