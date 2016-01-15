//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "ElFF/BBA05ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;

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
  gmp *= fMuP;
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
  gmn *= fMuN;
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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- load the BBA2005 fit coefficients
  fGep.a0 = fConfig->GetDoubleDef("Gep-a0", gc->GetDouble("BBA05-Gep-a0"));
  fGep.a1 = fConfig->GetDoubleDef("Gep-a1", gc->GetDouble("BBA05-Gep-a1"));
  fGep.a2 = fConfig->GetDoubleDef("Gep-a2", gc->GetDouble("BBA05-Gep-a2"));
  fGep.b1 = fConfig->GetDoubleDef("Gep-b1", gc->GetDouble("BBA05-Gep-b1"));
  fGep.b2 = fConfig->GetDoubleDef("Gep-b2", gc->GetDouble("BBA05-Gep-b2"));
  fGep.b3 = fConfig->GetDoubleDef("Gep-b3", gc->GetDouble("BBA05-Gep-b3"));
  fGep.b4 = fConfig->GetDoubleDef("Gep-b4", gc->GetDouble("BBA05-Gep-b4"));
  fGmp.a0 = fConfig->GetDoubleDef("Gmp-a0", gc->GetDouble("BBA05-Gmp-a0"));
  fGmp.a1 = fConfig->GetDoubleDef("Gmp-a1", gc->GetDouble("BBA05-Gmp-a1"));
  fGmp.a2 = fConfig->GetDoubleDef("Gmp-a2", gc->GetDouble("BBA05-Gmp-a2"));
  fGmp.b1 = fConfig->GetDoubleDef("Gmp-b1", gc->GetDouble("BBA05-Gmp-b1"));
  fGmp.b2 = fConfig->GetDoubleDef("Gmp-b2", gc->GetDouble("BBA05-Gmp-b2"));
  fGmp.b3 = fConfig->GetDoubleDef("Gmp-b3", gc->GetDouble("BBA05-Gmp-b3"));
  fGmp.b4 = fConfig->GetDoubleDef("Gmp-b4", gc->GetDouble("BBA05-Gmp-b4"));
  fGen.a0 = fConfig->GetDoubleDef("Gen-a0", gc->GetDouble("BBA05-Gen-a0"));
  fGen.a1 = fConfig->GetDoubleDef("Gen-a1", gc->GetDouble("BBA05-Gen-a1"));
  fGen.a2 = fConfig->GetDoubleDef("Gen-a2", gc->GetDouble("BBA05-Gen-a2"));
  fGen.b1 = fConfig->GetDoubleDef("Gen-b1", gc->GetDouble("BBA05-Gen-b1"));
  fGen.b2 = fConfig->GetDoubleDef("Gen-b2", gc->GetDouble("BBA05-Gen-b2"));
  fGen.b3 = fConfig->GetDoubleDef("Gen-b3", gc->GetDouble("BBA05-Gen-b3"));
  fGen.b4 = fConfig->GetDoubleDef("Gen-b4", gc->GetDouble("BBA05-Gen-b4"));
  fGmn.a0 = fConfig->GetDoubleDef("Gmn-a0", gc->GetDouble("BBA05-Gmn-a0"));
  fGmn.a1 = fConfig->GetDoubleDef("Gmn-a1", gc->GetDouble("BBA05-Gmn-a1"));
  fGmn.a2 = fConfig->GetDoubleDef("Gmn-a2", gc->GetDouble("BBA05-Gmn-a2"));
  fGmn.b1 = fConfig->GetDoubleDef("Gmn-b1", gc->GetDouble("BBA05-Gmn-b1"));
  fGmn.b2 = fConfig->GetDoubleDef("Gmn-b2", gc->GetDouble("BBA05-Gmn-b2"));
  fGmn.b3 = fConfig->GetDoubleDef("Gmn-b3", gc->GetDouble("BBA05-Gmn-b3"));
  fGmn.b4 = fConfig->GetDoubleDef("Gmn-b4", gc->GetDouble("BBA05-Gmn-b4"));

  //-- anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));
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
