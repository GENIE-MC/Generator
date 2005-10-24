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
double DipoleELFormFactorsModel::Gep(double q2) const
{
  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  double Mv2 = fConfig->GetDoubleDef("Mv2", kElMv2);

  // calculate and return GNE
  double ge = 1. / TMath::Power(1-q2/Mv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gep(q^2 = " << q2 << ") = " << ge;
  return ge;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gen(double /*q2*/) const
{
  return 0.;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmp(double q2) const
{
  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  double Mv2 = fConfig->GetDoubleDef("Mv2", kElMv2);

  // calculate & return Gm
  double gm = kMuP / TMath::Power(1-q2/Mv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmp(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gmn(double q2) const
{
  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  double Mv2 = fConfig->GetDoubleDef("Mv2", kElMv2);

  // calculate & return Gm
  double gm = kMuN / TMath::Power(1-q2/Mv2, 2);

  LOG("ELFormFactors", pDEBUG) << "Gmn(q^2 = " << q2 << ") = " << gm;
  return gm;
}
//____________________________________________________________________________

