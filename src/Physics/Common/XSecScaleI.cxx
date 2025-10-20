//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Physics/Common/XSecScaleI.h"

using namespace genie;

//_________________________________________________________________________
XSecScaleI::XSecScaleI( string name, string config /*"Default"*/ ) :
  Algorithm(name, config) 
{

}

//_________________________________________________________________________
XSecScaleI::~XSecScaleI()
{

}
//_________________________________________________________________________
void XSecScaleI::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig(); 
}
//____________________________________________________________________________
void XSecScaleI::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig(); 
}
//_________________________________________________________________________
