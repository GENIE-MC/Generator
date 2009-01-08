//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

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



