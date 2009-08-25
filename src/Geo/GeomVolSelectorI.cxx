//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

//#include "EVGDrivers/GeomVolSelectorI.h"
//rwh -- for now during test phase
#include "Geo/GeomVolSelectorI.h"

using namespace genie;
using namespace genie::geometry;

//____________________________________________________________________________
GeomVolSelectorI::GeomVolSelectorI() 
{
  fName = "no-name";
}

//____________________________________________________________________________
GeomVolSelectorI::GeomVolSelectorI(std::string name) 
{
  fName = name;
}

//___________________________________________________________________________
GeomVolSelectorI::~GeomVolSelectorI()
{

}
//___________________________________________________________________________
