//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionListGeneratorI.h"

using namespace genie;

//___________________________________________________________________________
InteractionListGeneratorI::InteractionListGeneratorI() :
Algorithm()
{

}
//___________________________________________________________________________
InteractionListGeneratorI::InteractionListGeneratorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
InteractionListGeneratorI::InteractionListGeneratorI(
                                                string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
InteractionListGeneratorI::~InteractionListGeneratorI()
{

}
//___________________________________________________________________________
