//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModelNC

\brief    Concrete implementation of the QELFormFactorsModelI :
          Form Factors for Quasi Elastic NC vN scattering according to
          Llewellyn-Smith model

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/ 
//____________________________________________________________________________

#include <iostream>

#include "Conventions/Constants.h"
#include "LlewellynSmith/LlewellynSmithModelNC.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LlewellynSmithModelNC::LlewellynSmithModelNC() :
LlewellynSmithModel()
{
  fName = "genie::LlewellynSmithModelNC";

  //FindConfig();
}
//____________________________________________________________________________
LlewellynSmithModelNC::LlewellynSmithModelNC(const char * param_set) :
LlewellynSmithModel(param_set)
{
  fName = "genie::LlewellynSmithModelNC";

  FindConfig();
}
//____________________________________________________________________________
LlewellynSmithModelNC::~LlewellynSmithModelNC()
{

}
//____________________________________________________________________________
double LlewellynSmithModelNC::F1V(const Interaction * interaction) const
{
  //-- calculate F1V-CC & Fp1
  
  double F1V_CC = LlewellynSmithModel::F1V(interaction);
  
  double Fp1    = LlewellynSmithModel::F1N(interaction);

  //-- calculate F1V-NC
  
  double F1V_NC = 0.5*F1V_CC - 2*kSin8w*kSin8w*Fp1;
  
  return F1V_NC;  
}
//____________________________________________________________________________
double LlewellynSmithModelNC::xiF2V(const Interaction * interaction) const
{
  //-- calculate xiF2V_CC and Fp2
  
  double xiF2V_CC = LlewellynSmithModel::xiF2V(interaction);

  double Fp2      = LlewellynSmithModel::munF2N(interaction) / kMuP;

  //-- calculate xiF2-NC
  
  double xiF2V_NC = 0.5*xiF2V_CC - 2*kSin8w*kSin8w*(kMuP-1)*Fp2;
  
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
  //-- get scattering parameters & auxiliary parameters

  const ScatteringParams & sc_params = interaction -> GetScatteringParams();

  double q2    = sc_params.q2();
  double Mnuc2 = pow(kNucleonMass, 2); // or init_state.TargetMass()? p/n
  double Mpi2  = pow(kPionMass,    2); 

  //-- calculate and return Fp
  
  double Fp_NC = 2 * Mnuc2 * ( this->FA(interaction) ) / ( Mpi2 - q2 );

  return Fp_NC;
}
//____________________________________________________________________________

