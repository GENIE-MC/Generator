//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithIsoFFCC.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithIsoFFCC::LwlynSmithIsoFFCC() :
LwlynSmithFF("genie::LwlynSmithIsoFFCC")
{

}
//____________________________________________________________________________
LwlynSmithIsoFFCC::LwlynSmithIsoFFCC(string config) :
LwlynSmithFF("genie::LwlynSmithIsoFFCC", config)
{

}
//____________________________________________________________________________
LwlynSmithIsoFFCC::~LwlynSmithIsoFFCC()
{

}
//____________________________________________________________________________
double LwlynSmithIsoFFCC::F1V(const Interaction * interaction) const
{
  if (fIsCC)
    return LwlynSmithFF::F1V(interaction);
    
  double F1p = this->F1P(interaction);
  double F1n = this->F1N(interaction);
  double _F1V = F1p + F1n;
  return _F1V;
}
//____________________________________________________________________________
double LwlynSmithIsoFFCC::xiF2V(const Interaction * interaction) const
{
  if (fIsCC)
    return LwlynSmithFF::xiF2V(interaction);
  
  double F2p = this->F2P(interaction);
  double F2n = this->F2N(interaction);
  double _xiF2V = F2p + F2n;
  return _xiF2V;
}
//____________________________________________________________________________
double LwlynSmithIsoFFCC::FA(const Interaction * interaction) const
{
  return LwlynSmithFF::FA(interaction);
}
//____________________________________________________________________________
double LwlynSmithIsoFFCC::Fp(const Interaction * interaction) const
{
  // get momentum transfer
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  // get struck nucleon mass & set pion mass
  PDGLibrary * pdglib = PDGLibrary::Instance();
  double M   = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;
  double M2  = M*M;
  double Mpi  = kPionMass;
  double Mpi2 = TMath::Power(Mpi, 2);

  // calculate FA
  double fa = this->FA(interaction);

  // calculate Fp
  double _Fp = 2. * M2 * fa/(Mpi2-q2);
  return _Fp;
}
//____________________________________________________________________________
double LwlynSmithIsoFFCC::tau(const Interaction * interaction) const
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
//____________________________________________________________________________
void LwlynSmithIsoFFCC::LoadConfig(void)
{
  LwlynSmithFF::LoadConfig();
  GetParamDef( "IsCCFormFactors", fIsCC, true);
}


