//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool - December 10, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Physics/MuonEnergyLoss/MuELossI.h"

using namespace genie;
using namespace genie::mueloss;

//___________________________________________________________________________
MuELossI::MuELossI() :
Algorithm()
{

}
//___________________________________________________________________________
MuELossI::MuELossI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
MuELossI::MuELossI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
MuELossI::~MuELossI()
{

}
//___________________________________________________________________________

