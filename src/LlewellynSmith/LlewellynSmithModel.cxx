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
LlewellynSmithModel::LlewellynSmithModel(const char * param_set) :
QELFormFactorsModelI(param_set)
{

}
//____________________________________________________________________________
LlewellynSmithModel::~LlewellynSmithModel()
{

}
//____________________________________________________________________________
double LlewellynSmithModel::F1V(const Interaction * interaction) const
{
  const ScatteringParams & sc_params = interaction -> GetScatteringParams();

  double q2     = sc_params.q2();
  double q2_4M2 = q2_4Mnucl2(interaction);

  double F1V = ( GVE(q2) - q2_4M2 * GVM(q2) ) / (1-q2_4M2);

  return F1V;
}
//____________________________________________________________________________
double LlewellynSmithModel::xiF2V(const Interaction * interaction) const
{
  const ScatteringParams & sc_params = interaction -> GetScatteringParams();

  double q2     = sc_params.q2();
  double q2_4M2 = q2_4Mnucl2(interaction);

  double xiF2V = ( GVM(q2) - GVE(q2) ) / (1-q2_4M2);

  return xiF2V;
}
//____________________________________________________________________________
double LlewellynSmithModel::FA(const Interaction * interaction) const
{
  //-- get Ma2 (axial mass squared) and FA(q2=0) from the configuration registry

  double Ma2 = 0, FA_0 = 0;

  fConfig->Get("Axial-Mass^2", Ma2);
  fConfig->Get("FA(q2=0)",     FA_0);

  //-- get scattering parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();

  double q2 = sc_params.q2();

  //-- calculate FA(q2)

  double FA = ( FA_0 / ( (1-q2/Ma2) * (1-q2/Ma2) ) );

  return FA;
}
//____________________________________________________________________________
double LlewellynSmithModel::Fp(const Interaction * interaction) const
{
  //-- get scattering & initial state parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();
  const InitialState & init_state    = interaction -> GetInitialState();

  double q2     = sc_params.q2();
  double Mnucl  = init_state.GetTarget().StruckNucleonMass();
  double Mnucl2 = pow(Mnucl, 2);
  
  //-- calculate Fp

  double Fp = ( 2 * Mnucl2 * FA(interaction) / (kPionMass*kPionMass-q2) );

  return Fp;
}
//____________________________________________________________________________
double LlewellynSmithModel::F1N(const Interaction * interaction) const
{
  //-- get auxiliary parameters

  double GNE    = this->GNE(interaction);
  double GNM    = this->GNM(interaction);
  double q2_4M2 = this->q2_4Mnucl2(interaction);

  //-- calculate & return F1N

  double F1N = (GNE - GNM*q2_4M2)/(1-q2_4M2);

  return F1N;
}
//____________________________________________________________________________
double LlewellynSmithModel::munF2N(const Interaction * interaction) const
{
  //-- get auxiliary parameters

  double GNE    = this->GNE(interaction);
  double GNM    = this->GNM(interaction);
  double q2_4M2 = this->q2_4Mnucl2(interaction);

  //-- calculate and return munF2N;

  double munF2N = (GNM-GNE) / (1-q2_4M2);

  return munF2N;
}
//____________________________________________________________________________
double LlewellynSmithModel::q2_4Mnucl2(
                                        const Interaction * interaction) const
{
  //-- get scattering & initial state parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();
  const InitialState & init_state    = interaction -> GetInitialState();

  double q2     = sc_params.q2();
  double Mnucl  = init_state.GetTarget().StruckNucleonMass();
  double Mnucl2 = pow(Mnucl, 2); 

  //-- calculate q^2 / (4*Mnuc^2)

  return q2/(4*Mnucl2);
}
//____________________________________________________________________________
double LlewellynSmithModel::GNE(const Interaction * interaction) const
{
  //-- get scattering & initial state parameters

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  double q2     = sc_params.q2();

  //-- get Mv2 (vector mass squared) from the configuration registry

  double Mv2 = 0;

  fConfig->Get("Vector-Mass^2", Mv2);

  //-- calculate and return GNE

  double GNE = this->GNE0(interaction) / pow(1-q2/Mv2, 2);

  return GNE;
}
//____________________________________________________________________________
double LlewellynSmithModel::GNE0(const Interaction * interaction) const
{
  //-- get initial state information

  const InitialState & init_state = interaction->GetInitialState();

  bool isP = init_state.GetTarget().IsProton();
  bool isN = init_state.GetTarget().IsNeutron();

  //-- return GNE0

  if      (isP) return  1;
  else if (isN) return  0;
  else          return -1; // complain here
}
//____________________________________________________________________________
double LlewellynSmithModel::GNM(const Interaction * interaction) const
{
  //-- get scattering parameters & auxiliary parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();

  double q2    = sc_params.q2();

  //-- get Mv2 (vector mass squared) from the configuration registry

  double Mv2 = 0;

  fConfig->Get("Vector-Mass^2", Mv2);

  //-- calculate & return GNM

  double GNM = ( this->GNM0(interaction) ) / pow(1-q2/Mv2, 2);

  return GNM;
}
//____________________________________________________________________________
double LlewellynSmithModel::GNM0(const Interaction * interaction) const
{
  //-- get initial state information

  const InitialState & init_state = interaction->GetInitialState();

  bool isP     = init_state.GetTarget().IsProton();
  bool isN     = init_state.GetTarget().IsNeutron();

  //-- get MuP and MuN (proton and neutron anomalous magnetic moments) from
  //   the configuration registry

  double MuP = 0, MuN = 0;

  fConfig->Get("Proton-Anomalous-Magnetic-Moment",  MuP);
  fConfig->Get("Neutron-Anomalous-Magnetic-Moment", MuN);

  //-- compute and return GNM0
  
  if      (isP) return MuP;
  else if (isN) return MuN;
  else          return -1;  // complain here
}
//____________________________________________________________________________
double LlewellynSmithModel::GVE(double q2) const
{
  //-- get Mv2 (vector mass squared) from the configuration registry

  double Mv2 = 0;

  fConfig->Get("Vector-Mass^2", Mv2);

  //-- calculate GVE

  double GVE = 1 / pow(1-q2/Mv2, 2);

  return GVE;
}
//____________________________________________________________________________
double LlewellynSmithModel::GVM(double q2) const
{
  //-- get Mv2 (vector mass squared), MuP and MuN (proton and neutron
  //   anomalous magnetic moments) from the configuration registry

  double Mv2 = 0, MuP = 0, MuN = 0;

  fConfig->Get("Vector-Mass^2",                     Mv2);
  fConfig->Get("Proton-Anomalous-Magnetic-Moment",  MuP);
  fConfig->Get("Neutron-Anomalous-Magnetic-Moment", MuN);

  //-- calculate GVM

  double GVM = (1+MuP-MuN) / pow(1-q2/Mv2, 2);

  return GVM;
}
//____________________________________________________________________________

