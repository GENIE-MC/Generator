//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModel

\brief    Abstract Base Class:
          implements the QELFormFactorsModelI interface but can not be
          instantiated.

          Its sole purpose of existence is to transmit common implementation
          (related to the Llewellyn-Smith model for QEL vN scattering) to its
          concrete subclasses: LlewellynSmithModelCC, LlewellynSmithModelNC.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/ELFormFactors.h"
#include "Base/ELFormFactorsModelI.h"
#include "Conventions/Constants.h"
#include "LlewellynSmith/LlewellynSmithModel.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LlewellynSmithModel::LlewellynSmithModel() :
QELFormFactorsModelI()
{

}
//____________________________________________________________________________
LlewellynSmithModel::LlewellynSmithModel(string name) :
QELFormFactorsModelI(name)
{

}
//____________________________________________________________________________
LlewellynSmithModel::LlewellynSmithModel(string name, string config) :
QELFormFactorsModelI(name, config)
{

}
//____________________________________________________________________________
LlewellynSmithModel::~LlewellynSmithModel()
{

}
//____________________________________________________________________________
double LlewellynSmithModel::F1V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double F1V = (gve - t*gvm) / (1-t);
  return F1V;
}
//____________________________________________________________________________
double LlewellynSmithModel::xiF2V(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gve = this->GVE(interaction);
  double gvm = this->GVM(interaction);

  double xiF2V = (gvm-gve) / (1-t);
  return xiF2V;
}
//____________________________________________________________________________
double LlewellynSmithModel::FA(const Interaction * interaction) const
{
  // get Ma2 (axial mass^2) and FA(q2=0) from the configuration registry or
  // set defaults
  double Ma2 = fConfig->GetDoubleDef("Ma2", kQelMa2);
  double FA0 = fConfig->GetDoubleDef("FA0", kQelFA0);

  // get scattering parameters
  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  // calculate FA(q2)
  double dn = TMath::Power(1.-q2/Ma2, 2);
  double FA = FA0/dn;
  return FA;
}
//____________________________________________________________________________
double LlewellynSmithModel::Fp(const Interaction * interaction) const
{
  // get momentum transfer
  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  // get struck nucleon mass & set pion mass
  const InitialState & init_state = interaction->GetInitialState();
  double MN   = init_state.GetTarget().StruckNucleonMass();
  double MN2  = TMath::Power(MN, 2);
  double Mpi  = kPionMass;
  double Mpi2 = TMath::Power(Mpi, 2);

  // calculate FA
  double fa = this->FA(interaction);

  // calculate Fp
  double Fp = 2. * MN2 * fa/(Mpi2-q2);
  return Fp;
}
//____________________________________________________________________________
double LlewellynSmithModel::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * MNucl^2)

  //-- get scattering & initial state parameters
  const ScatteringParams & scp = interaction->GetScatteringParams();
  const InitialState & init_state = interaction->GetInitialState();
  double q2     = scp.q2();
  double Mnucl  = init_state.GetTarget().StruckNucleonMass();
  double Mnucl2 = TMath::Power(Mnucl, 2);

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*Mnucl2);
}
//____________________________________________________________________________
double LlewellynSmithModel::GVE(const Interaction * interaction) const
{
  //-- get the momentum transfer
  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  //-- compute elastic form factors
  const ELFormFactorsModelI * elffmodel =
                 dynamic_cast<const ELFormFactorsModelI *>
                                          (this->SubAlg("el-form-factors"));
  ELFormFactors elff;
  elff.SetModel(elffmodel);
  elff.Calculate(q2);

  //-- compute GVE using CVC
  double gve = elff.Gep() - elff.Gen();
  return gve;
}
//____________________________________________________________________________
double LlewellynSmithModel::GVM(const Interaction * interaction) const
{
  //-- get the momentum transfer
  const ScatteringParams & scp = interaction->GetScatteringParams();
  double q2 = scp.q2();

  //-- compute elastic form factors
  const ELFormFactorsModelI * elffmodel =
                 dynamic_cast<const ELFormFactorsModelI *>
                                          (this->SubAlg("el-form-factors"));
  ELFormFactors elff;
  elff.SetModel(elffmodel);
  elff.Calculate(q2);

  //-- compute GVM using CVC
  double gvm = elff.Gmp() - elff.Gmn();
  return gvm;
}
//____________________________________________________________________________

