//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
