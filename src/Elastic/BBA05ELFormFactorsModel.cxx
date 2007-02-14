//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 19, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Elastic/BBA05Constants.h"
#include "Elastic/BBA05ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA05ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA05ELFormFactorsModel", config)
{

}
//____________________________________________________________________________
BBA05ELFormFactorsModel::~BBA05ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gep(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gep = this->BBA05Fit(t,fGep);
  return gep;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gmp = this->BBA05Fit(t,fGmp);
  return gmp;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gen = this->BBA05Fit(t,fGen);
  return gen;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gmn = this->BBA05Fit(t,fGmn);
  return gmn;
}
//____________________________________________________________________________
void BBA05ELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA05ELFormFactorsModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA05ELFormFactorsModel::LoadConfig(void)
{
  //-- load the BBA05 model parameters
  this->LoadBBA2005Params();
}
//____________________________________________________________________________
void BBA05ELFormFactorsModel::LoadBBA2005Params(void)
{
// Fill private data members holding BBA2003 model parameters.
// Use defaults for all configuration parameters not given in the config file.

  // BBA2004 fit coefficients

  fGep.a0 = fConfig->GetDoubleDef("Gep-a0", bba2005::kGep_a0);
  fGep.a1 = fConfig->GetDoubleDef("Gep-a1", bba2005::kGep_a1);
  fGep.a2 = fConfig->GetDoubleDef("Gep-a2", bba2005::kGep_a2);
  fGep.b1 = fConfig->GetDoubleDef("Gep-b1", bba2005::kGep_b1);
  fGep.b2 = fConfig->GetDoubleDef("Gep-b2", bba2005::kGep_b2);
  fGep.b3 = fConfig->GetDoubleDef("Gep-b3", bba2005::kGep_b3);
  fGep.b4 = fConfig->GetDoubleDef("Gep-b4", bba2005::kGep_b4);

  fGmp.a0 = fConfig->GetDoubleDef("Gmp-a0", bba2005::kGmp_a0);
  fGmp.a1 = fConfig->GetDoubleDef("Gmp-a1", bba2005::kGmp_a1);
  fGmp.a2 = fConfig->GetDoubleDef("Gmp-a2", bba2005::kGmp_a2);
  fGmp.b1 = fConfig->GetDoubleDef("Gmp-b1", bba2005::kGmp_b1);
  fGmp.b2 = fConfig->GetDoubleDef("Gmp-b2", bba2005::kGmp_b2);
  fGmp.b3 = fConfig->GetDoubleDef("Gmp-b3", bba2005::kGmp_b3);
  fGmp.b4 = fConfig->GetDoubleDef("Gmp-b4", bba2005::kGmp_b4);

  fGen.a0 = fConfig->GetDoubleDef("Gen-a0", bba2005::kGen_a0);
  fGen.a1 = fConfig->GetDoubleDef("Gen-a1", bba2005::kGen_a1);
  fGen.a2 = fConfig->GetDoubleDef("Gen-a2", bba2005::kGen_a2);
  fGen.b1 = fConfig->GetDoubleDef("Gen-b1", bba2005::kGen_b1);
  fGen.b2 = fConfig->GetDoubleDef("Gen-b2", bba2005::kGen_b2);
  fGen.b3 = fConfig->GetDoubleDef("Gen-b3", bba2005::kGen_b3);
  fGen.b4 = fConfig->GetDoubleDef("Gen-b4", bba2005::kGen_b4);

  fGmn.a0 = fConfig->GetDoubleDef("Gmn-a0", bba2005::kGmn_a0);
  fGmn.a1 = fConfig->GetDoubleDef("Gmn-a1", bba2005::kGmn_a1);
  fGmn.a2 = fConfig->GetDoubleDef("Gmn-a2", bba2005::kGmn_a2);
  fGmn.b1 = fConfig->GetDoubleDef("Gmn-b1", bba2005::kGmn_b1);
  fGmn.b2 = fConfig->GetDoubleDef("Gmn-b2", bba2005::kGmn_b2);
  fGmn.b3 = fConfig->GetDoubleDef("Gmn-b3", bba2005::kGmn_b3);
  fGmn.b4 = fConfig->GetDoubleDef("Gmn-b4", bba2005::kGmn_b4);
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::BBA05Fit(
                                      double t, const BBA2005Fit_t & fp) const
{
  double t2 = TMath::Power(t, 2);
  double t3 = TMath::Power(t, 3);
  double t4 = TMath::Power(t2,2);

  double Gn = (fp.a0) + (fp.a1*t) + (fp.a2*t2);
  double Gd = 1 + (fp.b1*t) + (fp.b2*t2) + (fp.b3*t3) + (fp.b4*t4);

  double G = Gn/Gd;
  return G;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::tau(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  const Target &     target = interaction->InitState().Tgt();

  double q2 = kine.q2(); // momentum transfer, <0
  double M2 = TMath::Power(target.HitNucMass(),2); // Mnucl^2

  double t = -q2/(4*M2); 
  return t;
}
//____________________________________________________________________________
