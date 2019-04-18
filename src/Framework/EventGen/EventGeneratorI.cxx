//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
