//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelNC

\brief    Concrete implementation of the QELFormFactorsModelI :
          Form Factors for Quasi Elastic NC vN scattering according to
          Llewellyn-Smith model.

\ref      E.A.Paschos and J.Y.Yu, hep-ph/0107261

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Base/ELFormFactors.h"
#include "Base/ELFormFactorsModelI.h"
#include "Conventions/Constants.h"
#include "LlewellynSmith/LlewellynSmithModelNC.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LlewellynSmithModelNC::LlewellynSmithModelNC() :
LlewellynSmithModel("genie::LlewellynSmithModelNC")
{

}
//____________________________________________________________________________
LlewellynSmithModelNC::LlewellynSmithModelNC(string config) :
LlewellynSmithModel("genie::LlewellynSmithModelNC", config)
{

}
//____________________________________________________________________________
LlewellynSmithModelNC::~LlewellynSmithModelNC()
{

}
//____________________________________________________________________________
double LlewellynSmithModelNC::F1V(const Interaction * interaction) const
{
  //-- calculate F1V-CC
  double F1V_CC = LlewellynSmithModel::F1V(interaction);

  //-- calculate F1p (see hep-ph/0107261)
  fELFF.Calculate(interaction);
  double t   = LlewellynSmithModel::tau(interaction);
  double F1p = fELFF.Gep() - t * fELFF.Gmp();

  //-- calculate F1V-NC
  double F1V_NC = 0.5*F1V_CC - 2*fSin28w*F1p;
  return F1V_NC;
}
//____________________________________________________________________________
double LlewellynSmithModelNC::xiF2V(const Interaction * interaction) const
{
  //-- calculate xiF2V_CC
  double xiF2V_CC = LlewellynSmithModel::xiF2V(interaction);

  //-- calculate F2p (see hep-ph/0107261)
  fELFF.Calculate(interaction);
  double F2p = (fELFF.Gmp() - fELFF.Gep()) / fMuP;

  //-- calculate xiF2-NC
  double xiF2V_NC = 0.5*xiF2V_CC - 2*fSin28w*(fMuP-1)*F2p;
  return xiF2V_NC;
}
//____________________________________________________________________________
double LlewellynSmithModelNC::FA(const Interaction * interaction) const
{
  //-- calculate FA_CC(q2)
  double FA_CC = LlewellynSmithModel::FA(interaction);

  //-- calculate & return FA_NC(q2)
  double FA_NC = 0.5 * FA_CC;
  return FA_NC;
}
//____________________________________________________________________________
double LlewellynSmithModelNC::Fp(const Interaction * interaction) const
{
  //-- get the momentum transfer
  const Kinematics & kine = interaction->GetKinematics();
  double q2 = kine.q2();

  //-- get struck nucleon mass & pion pass
  const InitialState & init_state = interaction->GetInitialState();
  double MN   = init_state.GetTarget().StruckNucleonMass();
  double MN2  = TMath::Power(MN,        2);
  double Mpi2 = TMath::Power(kPionMass, 2);

  //-- calculate FA
  double fa = this->FA(interaction);

  //-- calculate and return Fp
  double Fp_NC = 2*MN2*fa/(Mpi2-q2);
  return Fp_NC;
}
//____________________________________________________________________________

