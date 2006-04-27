//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - December 10, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "MuELoss/MuELossI.h"

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

