//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"

#include "Physics/Coherent/XSection/COHDeltaCurrent.h"

using namespace genie;


COHDeltaCurrent::COHDeltaCurrent() :
  COHFormFactorI("genie::COHDeltaCurrent")
{

}
//____________________________________________________________________________
COHDeltaCurrent::COHDeltaCurrent(string config) :
  COHFormFactorI("genie::COHDeltaCurrent", config)
{

}
//____________________________________________________________________________
COHDeltaCurrent::~COHDeltaCurrent()
{

}
//____________________________________________________________________________
GTrace R( const Interaction * i, 
	  const COHFormFactorI * ff ) const {
  
  GTrace t ;
  return t ;
}
//____________________________________________________________________________
void COHDeltaCurrent::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDeltaCurrent::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDeltaCurrent::LoadConfig(void)
{


}
//____________________________________________________________________________
