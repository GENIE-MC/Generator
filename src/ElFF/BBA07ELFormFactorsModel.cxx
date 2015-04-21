//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 31, 2009 - CA
   Was first added in v2.5.1.
 @ Sep 19, 2009 - CA
   Moved into the ElFF package from its previous location               

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "ElFF/BBA07ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;
using namespace genie::constants;

// Some model parameters - hardcoded for the time-being
//
double BBA07ELFormFactorsModel::fsParGalsterFactor [3] = { 1.7, 3.3, kNeutronMass };
double BBA07ELFormFactorsModel::fsParGepKelly      [5] = { -0.24,   10.98, 12.82, 21.97, kProtonMass };
double BBA07ELFormFactorsModel::fsParGmpKelly      [5] = {  0.1717, 11.26, 19.32, 8.33,  kProtonMass };
double BBA07ELFormFactorsModel::fsParGepLagrange   [8] = { 1., 0.9927, 0.9898, 0.9975, 0.9812, 0.9340, 1.0000, kProtonMass  };
double BBA07ELFormFactorsModel::fsParGmpLagrange   [8] = { 1., 1.0011, 0.9992, 0.9974, 1.0010, 1.0003, 1.0000, kProtonMass  };
double BBA07ELFormFactorsModel::fsParGmnLagrange_25[8] = { 1., 0.9958, 0.9877, 1.0193, 1.0350, 0.9164, 0.7300, kNeutronMass };
double BBA07ELFormFactorsModel::fsParGmnLagrange_43[8] = { 1., 0.9959, 0.9851, 1.0187, 1.0307, 0.9080, 0.9557, kNeutronMass };
double BBA07ELFormFactorsModel::fsParGenLagrange_25[8] = { 1., 1.1011, 1.1392, 1.0203, 1.1093, 1.5429, 0.9706, kNeutronMass };
double BBA07ELFormFactorsModel::fsParGenLagrange_43[8] = { 1., 1.1019, 1.1387, 1.0234, 1.1046, 1.5395, 1.2708, kNeutronMass };

//____________________________________________________________________________
BBA07ELFormFactorsModel::BBA07ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA07ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA07ELFormFactorsModel::BBA07ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA07ELFormFactorsModel", config)
{

}
//____________________________________________________________________________
BBA07ELFormFactorsModel::~BBA07ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gep(const Interaction * interaction) const
{
  double lagr  = this->Lagrange(interaction, fsParGepLagrange);
  double kelly = this->Kelly(interaction, fsParGepKelly);
  double gep   = lagr*kelly;

  return gep;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  double lagr  = this->Lagrange(interaction, fsParGmpLagrange);
  double kelly = this->Kelly(interaction, fsParGmpKelly);
  double gmp   = lagr*kelly; 

  return gmp;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double lagr  = this->Lagrange(interaction, fsParGenLagrange_43);
  double galst = this->GalsterFactor(interaction, fsParGalsterFactor);
  double gep   = this->Gep(interaction);
  double gen   = lagr*galst*gep;

  return gen;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  double lagr = this->Lagrange(interaction, fsParGmnLagrange_43);
  double gmp  = this->Gmp(interaction);
  double gmn  = lagr*gmp;

  return gmn;
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));

  // other parameters hardcoded for the time-being
  // ...

}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Tau(const Interaction * interaction) const
{
// tau = Q2 / (4*M2)

  const Kinematics & kine   = interaction->Kine();
  const Target &     target = interaction->InitState().Tgt();

  double q2 = kine.q2(); // momentum transfer, <0
  double M2 = TMath::Power(target.HitNucMass(),2); // Mnucl^2

  double t = -q2/(4*M2); 
  return t;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Xi(const Interaction * interaction) const
{
// Nachtman variable xi

  double tau = this->Tau(interaction);
  if(tau!=0.) {
    return 2./(1.+TMath::Sqrt(1.+1./tau));
  } 
  return 0;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Lagrange(
          const Interaction * interaction, double* par) const
{
// Lagrange parameterization

  static const int N = 7;
  static const double nodes[N] = {0., 1./6., 2./6., 3./6., 4./6., 5./6., 1.};

  double xi = this->Xi(interaction);

  double sum = 0.;
  for (int i = 0; i<N; ++i) {
     double part = par[i];
     for (int j = 0; j<N; ++j) {
            if ( i != j ) {
                part *= (xi - nodes[j])/(nodes[i]-nodes[j]);
            }
        }
        sum += part;
  }
  return sum;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Kelly(
          const Interaction * interaction, double* par) const
{
// Kelly parameterization

  double tau = this->Tau(interaction);

  double numerator   = 1.;
  double denominator = 1.;
  numerator   += par[0]*tau;
  denominator += par[1]*tau;
  denominator += par[2]*tau*tau;
  denominator += par[3]*tau*tau*tau;

  double kelly = (denominator!=0.) ? numerator/denominator : 0.;
  return kelly;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::GalsterFactor(
          const Interaction * interaction, double* par) const
{
// "Galster factor" a*tau/(1+b*tau) 

  double tau = this->Tau(interaction);
  double gf  = par[0]*tau/(1. + par[1]*tau);
  return gf;
}
//____________________________________________________________________________

