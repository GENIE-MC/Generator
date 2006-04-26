//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "LlewellynSmith/LlewellynSmithModelCC.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LlewellynSmithModelCC::LlewellynSmithModelCC() :
LlewellynSmithModel("genie::LlewellynSmithModelCC")
{

}
//____________________________________________________________________________
LlewellynSmithModelCC::LlewellynSmithModelCC(string config) :
LlewellynSmithModel("genie::LlewellynSmithModelCC", config)
{

}
//____________________________________________________________________________
LlewellynSmithModelCC::~LlewellynSmithModelCC()
{

}
//____________________________________________________________________________
double LlewellynSmithModelCC::F1V(const Interaction * interaction) const
{
  return LlewellynSmithModel::F1V(interaction);
}
//____________________________________________________________________________
double LlewellynSmithModelCC::xiF2V(const Interaction * interaction) const
{
  return LlewellynSmithModel::xiF2V(interaction);
}
//____________________________________________________________________________
double LlewellynSmithModelCC::FA(const Interaction * interaction) const
{
  return LlewellynSmithModel::FA(interaction);
}
//____________________________________________________________________________
double LlewellynSmithModelCC::Fp(const Interaction * interaction) const
{
  return LlewellynSmithModel::Fp(interaction);
}
//____________________________________________________________________________



