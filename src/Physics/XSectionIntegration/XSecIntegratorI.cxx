//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;

//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
XSecIntegratorI::~XSecIntegratorI()
{

}
//___________________________________________________________________________
