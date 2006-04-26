//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Fragmentation/HadronizationModelI.h"

using namespace genie;

//____________________________________________________________________________
HadronizationModelI::HadronizationModelI() :
Algorithm()
{

}
//____________________________________________________________________________
HadronizationModelI::HadronizationModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
HadronizationModelI::HadronizationModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
HadronizationModelI::~HadronizationModelI() 
{ 

}
//____________________________________________________________________________



