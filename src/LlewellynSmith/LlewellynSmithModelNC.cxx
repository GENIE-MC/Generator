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

  //-- compute elastic form factors
  const ELFormFactorsModelI * elffmodel =
                 dynamic_cast<const ELFormFactorsModelI *>
                                          (this->SubAlg("el-form-factors"));
  ELFormFactors elff;
  elff.SetModel(elffmodel);

  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  elff.Calculate(q2);

  //-- calculate F1p (see hep-ph/0107261)
  double t   = LlewellynSmithModel::tau(interaction);
  double F1p = elff.Gep() - t * elff.Gmp();

  //-- calculate F1V-NC

  double w2 = kSin8w_2; // sin^2(weinberg-angle)
  double F1V_NC = 0.5*F1V_CC - 2*w2*F1p;
  return F1V_NC;
}
//____________________________________________________________________________
double LlewellynSmithModelNC::xiF2V(const Interaction * interaction) const
{
  //-- calculate xiF2V_CC
  double xiF2V_CC = LlewellynSmithModel::xiF2V(interaction);

  //-- compute elastic form factors
  const ELFormFactorsModelI * elffmodel =
                 dynamic_cast<const ELFormFactorsModelI *>
                                          (this->SubAlg("el-form-factors"));
  ELFormFactors elff;
  elff.SetModel(elffmodel);

  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  elff.Calculate(q2);

  //-- calculate F2p (see hep-ph/0107261)
  double F2p = (elff.Gmp() - elff.Gep()) / kMuP;

  //-- calculate xiF2-NC
  double w2 = kSin8w_2; // sin^2(weinberg-angle)
  double xiF2V_NC = 0.5*xiF2V_CC - 2*w2*(kMuP-1)*F2p;
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
  const ScatteringParams & sc_params = interaction->GetScatteringParams();
  double q2 = sc_params.q2();

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

