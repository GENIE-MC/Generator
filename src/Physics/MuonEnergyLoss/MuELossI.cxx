//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - December 10, 2003

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

