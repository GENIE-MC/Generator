//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
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
