//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Framework/EventGen/EventGeneratorI.h"

using namespace genie;

//___________________________________________________________________________
EventGeneratorI::EventGeneratorI() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
EventGeneratorI::EventGeneratorI(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
EventGeneratorI::EventGeneratorI(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
EventGeneratorI::~EventGeneratorI()
{

}
//___________________________________________________________________________
