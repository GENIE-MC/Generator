//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - June 15, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
