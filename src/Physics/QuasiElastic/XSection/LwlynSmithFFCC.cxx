//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Renamed LlewellynSmithModelCC -> LwlynSmithFFCC

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



