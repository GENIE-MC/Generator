//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "ReinSehgal/BergerSehgalRESPXSec2014.h"

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
