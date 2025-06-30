//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <komninos-john.plows \at physics.ox.ac.uk>
 University of Oxford
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/HNLChannelCalculatorI.h"

using namespace genie;
using namespace genie::hnl;

//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI() : Algorithm()
{

}
//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI(string name) : Algorithm(name)
{

}
//____________________________________________________________________________
ChannelCalculatorI::ChannelCalculatorI(string name, string config) : 
  Algorithm(name, config)
{

}
//____________________________________________________________________________
