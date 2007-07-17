//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Base/XSecIntegratorI.h"

using namespace genie;

//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
XSecIntegratorI::XSecIntegratorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
XSecIntegratorI::~XSecIntegratorI()
{

}
//___________________________________________________________________________
