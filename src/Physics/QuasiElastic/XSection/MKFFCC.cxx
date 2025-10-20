//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/MKFFCC.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MKFFCC::MKFFCC() :
LwlynSmithFF("genie::MKFFCC")
{

}
//____________________________________________________________________________
MKFFCC::MKFFCC(string config) :
LwlynSmithFF("genie::MKFFCC", config)
{

}
//____________________________________________________________________________
MKFFCC::~MKFFCC()
{

}
//____________________________________________________________________________
double MKFFCC::F1V(const Interaction * interaction) const
{
  return LwlynSmithFF::F1V(interaction);
}
//____________________________________________________________________________
double MKFFCC::xiF2V(const Interaction * interaction) const
{
  return LwlynSmithFF::xiF2V(interaction);
}
//____________________________________________________________________________
double MKFFCC::FA(const Interaction * interaction) const
{
  return LwlynSmithFF::FA(interaction);
}
//____________________________________________________________________________
double MKFFCC::Fp(const Interaction * interaction) const
{
  return LwlynSmithFF::Fp(interaction);
}
//____________________________________________________________________________
double MKFFCC::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * Misoscalar^2)

  //-- get kinematics & initial state parameters
  const Kinematics &   kinematics = interaction->Kine();
  //const InitialState & init_state = interaction->InitState();
  double q2     = kinematics.q2();
 
  PDGLibrary * pdglib = PDGLibrary::Instance();
  double M = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*M*M);
}



