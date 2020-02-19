//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/Hadronization/FragmentationFunctionI.h"

using namespace genie;

//___________________________________________________________________________
FragmentationFunctionI::FragmentationFunctionI() :
Algorithm()
{

}
//___________________________________________________________________________
FragmentationFunctionI::FragmentationFunctionI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
FragmentationFunctionI::FragmentationFunctionI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
FragmentationFunctionI::~FragmentationFunctionI()
{

}
//___________________________________________________________________________
