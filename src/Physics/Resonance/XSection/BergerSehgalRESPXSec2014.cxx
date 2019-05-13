//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/BergerSehgalRESPXSec2014.h"

using namespace genie;

//____________________________________________________________________________
BergerSehgalRESPXSec2014::BergerSehgalRESPXSec2014() :
BSKLNBaseRESPXSec2014("genie::BergerSehgalRESPXSec2014")
{
  this->fKLN = false;
  this->fBRS = true;
}
//____________________________________________________________________________
BergerSehgalRESPXSec2014::BergerSehgalRESPXSec2014(string config) :
BSKLNBaseRESPXSec2014("genie::BergerSehgalRESPXSec2014", config)
{
  this->fKLN = false;
  this->fBRS = true;
}
//____________________________________________________________________________
BergerSehgalRESPXSec2014::~BergerSehgalRESPXSec2014()
{

}
//____________________________________________________________________________
