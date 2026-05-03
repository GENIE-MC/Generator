//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
*/
//____________________________________________________________________________

#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithFFEM.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFFEM::LwlynSmithFFEM() :
LwlynSmithFF("genie::LwlynSmithFFEM")
{

}
//____________________________________________________________________________
LwlynSmithFFEM::LwlynSmithFFEM(string config) :
LwlynSmithFF("genie::LwlynSmithFFEM", config)
{

}
//____________________________________________________________________________
LwlynSmithFFEM::~LwlynSmithFFEM()
{

}
//____________________________________________________________________________
double LwlynSmithFFEM::F1V(const Interaction* interaction) const
{
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  if ( pdg::IsProton(hit_nuc_pdg) ) return LwlynSmithFF::F1P(interaction);
  else return LwlynSmithFF::F1N(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFEM::xiF2V(const Interaction* interaction) const
{
  int hit_nuc_pdg = interaction->InitState().Tgt().HitNucPdg();
  if ( pdg::IsProton(hit_nuc_pdg) ) return LwlynSmithFF::F2P(interaction);
  else return LwlynSmithFF::F2N(interaction);
}
