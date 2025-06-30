//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
*/
//____________________________________________________________________________

#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithFFCC.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFFCC::LwlynSmithFFCC() :
LwlynSmithFF("genie::LwlynSmithFFCC")
{

}
//____________________________________________________________________________
LwlynSmithFFCC::LwlynSmithFFCC(string config) :
LwlynSmithFF("genie::LwlynSmithFFCC", config)
{

}
//____________________________________________________________________________
LwlynSmithFFCC::~LwlynSmithFFCC()
{

}
//____________________________________________________________________________
double LwlynSmithFFCC::F1V(const Interaction * interaction) const
{
  return LwlynSmithFF::F1V(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFCC::xiF2V(const Interaction * interaction) const
{
  return LwlynSmithFF::xiF2V(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFCC::FA(const Interaction * interaction) const
{
  return LwlynSmithFF::FA(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFCC::Fp(const Interaction * interaction) const
{
  return LwlynSmithFF::Fp(interaction);
}
//____________________________________________________________________________
