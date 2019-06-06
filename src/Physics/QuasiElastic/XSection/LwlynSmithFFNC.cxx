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
   Renamed LlewellynSmithModelNC -> LwlynSmithFFNC

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/QuasiElastic/XSection/ELFormFactors.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/LwlynSmithFFNC.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFFNC::LwlynSmithFFNC() :
LwlynSmithFF("genie::LwlynSmithFFNC")
{

}
//____________________________________________________________________________
LwlynSmithFFNC::LwlynSmithFFNC(string config) :
LwlynSmithFF("genie::LwlynSmithFFNC", config)
{

}
//____________________________________________________________________________
LwlynSmithFFNC::~LwlynSmithFFNC()
{

}
//____________________________________________________________________________
double LwlynSmithFFNC::F1V(const Interaction * interaction) const
{
  //-- calculate F1V-CC
  double F1V_CC = LwlynSmithFF::F1V(interaction);

  //-- calculate F1p (see hep-ph/0107261)
  fELFF.Calculate(interaction);
  double t   = LwlynSmithFF::tau(interaction);
  double F1p = fELFF.Gep() - t * fELFF.Gmp();

  //-- calculate F1V-NC
  double F1V_NC = 0.5*F1V_CC - 2*fSin28w*F1p;
  return F1V_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::xiF2V(const Interaction * interaction) const
{
  //-- calculate xiF2V_CC
  double xiF2V_CC = LwlynSmithFF::xiF2V(interaction);

  //-- calculate F2p (see hep-ph/0107261)
  fELFF.Calculate(interaction);
  double F2p = (fELFF.Gmp() - fELFF.Gep()) / fMuP;

  //-- calculate xiF2-NC
  double xiF2V_NC = 0.5*xiF2V_CC - 2*fSin28w*(fMuP-1)*F2p;
  return xiF2V_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::FA(const Interaction * interaction) const
{
  //-- calculate FA_CC(q2)
  double FA_CC = LwlynSmithFF::FA(interaction);

  //-- calculate & return FA_NC(q2)
  double FA_NC = 0.5 * FA_CC;
  return FA_NC;
}
//____________________________________________________________________________
double LwlynSmithFFNC::Fp(const Interaction * interaction) const
{
  //-- get the momentum transfer
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  //-- get struck nucleon mass & pion pass
  const InitialState & init_state = interaction->InitState();
  double MN   = init_state.Tgt().HitNucMass();
  double MN2  = TMath::Power(MN,        2);
  double Mpi2 = TMath::Power(kPionMass, 2);

  //-- calculate FA
  double fa = this->FA(interaction);

  //-- calculate and return Fp
  double Fp_NC = 2*MN2*fa/(Mpi2-q2);
  return Fp_NC;
}
//____________________________________________________________________________

