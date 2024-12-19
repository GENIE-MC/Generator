//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/EventGen/EventRecordVisitorI.h"

using namespace genie;

//___________________________________________________________________________
EventRecordVisitorI::EventRecordVisitorI() :
Algorithm()
{

}
//___________________________________________________________________________
EventRecordVisitorI::EventRecordVisitorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
EventRecordVisitorI::EventRecordVisitorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
EventRecordVisitorI::~EventRecordVisitorI()
{

}
//___________________________________________________________________________
