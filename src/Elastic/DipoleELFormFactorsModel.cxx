//____________________________________________________________________________
/*!

\class    genie::DipoleELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes dipole elastic form factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 19, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Elastic/DipoleELFormFactorsModel.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DipoleELFormFactorsModel::DipoleELFormFactorsModel() :
ELFormFactorsModelI("genie::DipoleELFormFactorsModel")
{

}
//____________________________________________________________________________
DipoleELFormFactorsModel::DipoleELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::DipoleELFormFactorsModel", config)
{

}
//____________________________________________________________________________
DipoleELFormFactorsModel::~DipoleELFormFactorsModel()
{

}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gep(const Interaction * interaction) const
{
  // calculate and return GNE
  double q2 = interaction->GetKinematics().q2();
  double ge = 1. / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gep(q^2 = " << q2 << ") = " << ge;
  return ge;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gen(const Interaction * /*in*/) const
{
  return 0.;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  // calculate & return Gm
  double q2 = interaction->GetKinematics().q2();
  double gm = kMuP / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmp(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  // calculate & return Gm
  double q2 = interaction->GetKinematics().q2();
  double gm = kMuN / TMath::Power(1-q2/fMv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmn(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void DipoleELFormFactorsModel::LoadConfig(void)
{
  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  fMv2 = fConfig->GetDoubleDef("Mv2", kElMv2);
}
//____________________________________________________________________________

