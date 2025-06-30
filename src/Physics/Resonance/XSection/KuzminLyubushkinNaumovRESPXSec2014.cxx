//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steve Dytman
 University of Pittsburgh

 Jarek Nowak
 University of Lancaster

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/KuzminLyubushkinNaumovRESPXSec2014.h"

using namespace genie;

//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::KuzminLyubushkinNaumovRESPXSec2014() :
BSKLNBaseRESPXSec2014("genie::KuzminLyubushkinNaumovRESPXSec2014")
{
  this->fKLN = true;
  this->fBRS = false;
}
//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::KuzminLyubushkinNaumovRESPXSec2014(string config) :
BSKLNBaseRESPXSec2014("genie::KuzminLyubushkinNaumovRESPXSec2014", config)
{
  this->fKLN = true;
  this->fBRS = false;
}
//____________________________________________________________________________
KuzminLyubushkinNaumovRESPXSec2014::~KuzminLyubushkinNaumovRESPXSec2014()
{

}
//____________________________________________________________________________
