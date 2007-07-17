//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Numerical/MinimizerI.h"

using namespace genie;

//____________________________________________________________________________
MinimizerI::MinimizerI() :
Algorithm()
{

}
//____________________________________________________________________________
MinimizerI::MinimizerI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
MinimizerI::MinimizerI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
MinimizerI::~MinimizerI()
{

}
//____________________________________________________________________________
