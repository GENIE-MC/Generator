//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "EVGCore/EventRecordVisitorI.h"

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
