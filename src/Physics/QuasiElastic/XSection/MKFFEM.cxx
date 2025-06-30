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

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/MKFFEM.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MKFFEM::MKFFEM() :
QELFormFactorsModelI("genie::MKFFEM")
{

}
//____________________________________________________________________________
MKFFEM::MKFFEM(string config) :
QELFormFactorsModelI("genie::MKFFEM", config)
{

}
//____________________________________________________________________________
MKFFEM::~MKFFEM()
{

}
//____________________________________________________________________________
double MKFFEM::F1P(const Interaction * interaction) const
{ 
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gep() - t * fELFF.Gmp());
}
//____________________________________________________________________________
double MKFFEM::F2P(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmp() - fELFF.Gep());
}
//____________________________________________________________________________
double MKFFEM::F1N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gen() - t * fELFF.Gmn());
}
//____________________________________________________________________________
double MKFFEM::F2N(const Interaction * interaction) const
{
  fELFF.Calculate(interaction);
  double t   = this->tau(interaction);
  double T   = 1 / (1 - t);
  return T * (fELFF.Gmn() - fELFF.Gen());
}
//____________________________________________________________________________
double MKFFEM::F1V(const Interaction * interaction) const
{
  double F1p = this->F1P(interaction);
  double F1n = this->F1N(interaction);

  double _F1V = F1p + F1n;
  return _F1V;
}
//____________________________________________________________________________
double MKFFEM::xiF2V(const Interaction * interaction) const
{
  double F2p = this->F2P(interaction);
  double F2n = this->F2N(interaction);

  double _xiF2V = F2p + F2n;
  return _xiF2V;
}
//____________________________________________________________________________
double MKFFEM::FA(const Interaction * /*interaction*/ ) const
{
  return 0.;
}
//____________________________________________________________________________
double MKFFEM::Fp(const Interaction * /*interaction*/ ) const
{
  return 0.;
}
//____________________________________________________________________________
void MKFFEM::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKFFEM::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MKFFEM::LoadConfig(void)
{
// Load configuration data from its configuration Registry (or global defaults)
// to private data members
  fElFFModel =
    dynamic_cast<const ELFormFactorsModelI *> (this->SubAlg("ElasticFormFactorsModel"));
  assert(fElFFModel);
  fELFF.SetModel(fElFFModel);
  

}
//____________________________________________________________________________
double MKFFEM::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * MNucl^2)

  //-- get kinematics & initial state parameters
  const Kinematics &   kinematics = interaction->Kine();
  //  const InitialState & init_state = interaction->InitState();
  double q2     = kinematics.q2();
  
  PDGLibrary * pdglib = PDGLibrary::Instance();
  double M = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*M*M);
}
//____________________________________________________________________________
