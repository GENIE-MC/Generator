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
  // get scattering parameters
  const Kinematics & kine = interaction->GetKinematics();
  double q2 = kine.q2();

  // calculate FA(q2)
  double dn = TMath::Power(1.-q2/fMa2, 2);
  double FA = fFA0/dn;
  return FA;
}
//____________________________________________________________________________
double LlewellynSmithModel::Fp(const Interaction * interaction) const
{
  // get momentum transfer
  const Kinematics & kine = interaction->GetKinematics();
  double q2 = kine.q2();

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
void LlewellynSmithModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->LoadConfigData();
}
//____________________________________________________________________________
void LlewellynSmithModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
  this->LoadConfigData();
}
//____________________________________________________________________________
void LlewellynSmithModel::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fElFFModel = 0;

  //-- load elastic form factors model
  fElFFModel = dynamic_cast<const ELFormFactorsModelI *>
                                          (this->SubAlg("el-form-factors"));
  assert(fElFFModel);
}
//____________________________________________________________________________
void LlewellynSmithModel::LoadConfigData(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.
// Sets defaults for configuration options that were not specified

  fMa2 = fConfig->GetDoubleDef("fMa2", kQelMa2); // fMa2 (axial mass^2)
  fFA0 = fConfig->GetDoubleDef("fFA0", kQelFA0); // FA(q2=0)

  assert(fMa2>0);

  fELFF.SetModel(fElFFModel);  
}
//____________________________________________________________________________
double LlewellynSmithModel::tau(const Interaction * interaction) const
{
// computes q^2 / (4 * MNucl^2)

  //-- get kinematics & initial state parameters
  const Kinematics &   kinematics = interaction->GetKinematics();
  const InitialState & init_state = interaction->GetInitialState();
  double q2     = kinematics.q2();
  double Mnucl  = init_state.GetTarget().StruckNucleonMass();
  double Mnucl2 = TMath::Power(Mnucl, 2);

  //-- calculate q^2 / (4*Mnuc^2)
  return q2/(4*Mnucl2);
}
//____________________________________________________________________________
double LlewellynSmithModel::GVE(const Interaction * interaction) const
{
  //-- compute GVE using CVC

  fELFF.Calculate(interaction);
  double gve = fELFF.Gep() - fELFF.Gen();
  return gve;
}
//____________________________________________________________________________
double LlewellynSmithModel::GVM(const Interaction * interaction) const
{
  //-- compute GVM using CVC

  fELFF.Calculate(interaction);
  double gvm = fELFF.Gmp() - fELFF.Gmn();
  return gvm;
}
//____________________________________________________________________________

