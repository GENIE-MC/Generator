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
ELFormFactorsModelI()
{

}
//____________________________________________________________________________
DipoleELFormFactorsModel::DipoleELFormFactorsModel(const char * param_set) :
ELFormFactorsModelI(param_set)
{
  fName = "genie::DipoleELFormFactorsModel";

  this->FindConfig();
}
//____________________________________________________________________________
DipoleELFormFactorsModel::~DipoleELFormFactorsModel()
{

}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Ge(const Interaction * interaction) const
{
  // get initial state & scattering parameters
  const InitialState & init_state = interaction->GetInitialState();
  const Target & target = init_state.GetTarget();
  const ScatteringParams & scp = interaction -> GetScatteringParams();

  // get Ge(0) for the struck nucleon
  double ge0;
  if      ( target.IsProton()  ) ge0 = 1.;
  else if ( target.IsNeutron() ) ge0 = 0.;
  else  {
    LOG("ELFormFactors", pERROR) << "Undefined struck nucleon!";
    return -99999;
  }

  // get the momentum transfer
  double q2 = scp.q2();

  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  double Mv2 = fConfig->GetDoubleDef("Mv2", kElMv2);

  // calculate and return GNE
  double ge = ge0 / TMath::Power(1-q2/Mv2, 2);

  return ge;
}
//____________________________________________________________________________
double DipoleELFormFactorsModel::Gm(const Interaction * interaction) const
{
  // get initial state & scattering parameters
  const InitialState & init_state = interaction->GetInitialState();
  const Target & target = init_state.GetTarget();
  const ScatteringParams & scp = interaction->GetScatteringParams();

  // get the anomalous magnetic moment of the struck nucleon
  double MagnMom;
  if      ( target.IsProton()  ) MagnMom = kMuP;
  else if ( target.IsNeutron() ) MagnMom = kMuN;
  else  {
    LOG("ELFormFactors", pERROR) << "Undefined struck nucleon!";
    return -99999;
  }

  // get the momentum transfer
  double q2 = scp.q2();

  // get Mv2 (vector mass squared) from the configuration registry if it
  // exists, otherwise use the default value
  double Mv2 = fConfig->GetDoubleDef("Mv2", kElMv2);

  // calculate & return Gm
  double Gm = MagnMom / TMath::Power(1-q2/Mv2, 2);

  return Gm;
}
//____________________________________________________________________________

