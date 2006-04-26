//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 19, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Elastic/BBA03Constants.h"
#include "Elastic/BBA03ELFormFactorsModel.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA03ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA03ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA03ELFormFactorsModel::~BBA03ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gep(const Interaction * interaction) const
{
  double gep = 0;
  double q2  = interaction->GetKinematics().q2();

  if( TMath::Abs(q2) > fQ2Max ) {
     double gepmx = this->BBA03Fit(-fQ2Max, 1.,   fGep);
     double gmpmx = this->BBA03Fit(-fQ2Max, fMuP, fGmp);
     double gmp   = this->BBA03Fit(q2, fMuP, fGmp);
     gep = gmp * (gepmx/gmpmx);
  } else {
     gep = this->BBA03Fit(q2, 1., fGep);
  }
  return gep;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  double q2  = interaction->GetKinematics().q2();
  double gmp = this->BBA03Fit(q2, fMuP, fGmp);
  return gmp;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double q2  = interaction->GetKinematics().q2();

  const Target & tgt = interaction->GetInitialState().GetTarget();

  double M   = tgt.StruckNucleonMass();      // Mnucl
  double M2  = TMath::Power(M,2);            // Mnucl^2
  double t   = -q2/(4*M2);                   // q2<0
  double a   = fGenA;                        // Krutov et al. parameter a
  double b   = fGenB;                        // Krutov et al. parameter b
  double mv2 = fMv2;                         // elastic vector mass^2
  double GD  = 1./TMath::Power(1-q2/mv2,2.); // dipole form factor

  double gen = -1. * fMuN * a * t * GD / (1 + b*t);
  return gen;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  double q2  = interaction->GetKinematics().q2();
  double gmn = this->BBA03Fit(q2, fMuN, fGmn);
  return gmn;
}
//____________________________________________________________________________
void BBA03ELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA03ELFormFactorsModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA03ELFormFactorsModel::LoadConfig(void)
{
  //-- load BBA03 model parameters

  this->LoadBBA2003Params();

  //-- get config options from the configuration registry or set defaults
  //   from the global parameter list

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // vector mass
  fMv  = fConfig->GetDoubleDef("Mv", gc->GetDouble("EL-Mv"));
  fMv2 = TMath::Power(fMv,2);

  // anomalous magnetic moments
  fMuP = fConfig->GetDoubleDef("MuP", gc->GetDouble("AnomMagnMoment-P"));
  fMuN = fConfig->GetDoubleDef("MuN", gc->GetDouble("AnomMagnMoment-N"));
}
//____________________________________________________________________________
void BBA03ELFormFactorsModel::LoadBBA2003Params(void)
{
// Fill private data members holding BBA2003 model parameters.
// Use defaults for all configuration parameters not given in the config file.

  // BBA2003 fit coefficients
  fGep.a2  = fConfig->GetDoubleDef("Gep-a2",  bba2003::kGep_a2 );
  fGep.a4  = fConfig->GetDoubleDef("Gep-a4",  bba2003::kGep_a4 );
  fGep.a6  = fConfig->GetDoubleDef("Gep-a6",  bba2003::kGep_a6 );
  fGep.a8  = fConfig->GetDoubleDef("Gep-a8",  bba2003::kGep_a8 );
  fGep.a10 = fConfig->GetDoubleDef("Gep-a10", bba2003::kGep_a10);
  fGep.a12 = fConfig->GetDoubleDef("Gep-a12", bba2003::kGep_a12);

  fGmp.a2  = fConfig->GetDoubleDef("Gmp-a2",  bba2003::kGmp_a2 );
  fGmp.a4  = fConfig->GetDoubleDef("Gmp-a4",  bba2003::kGmp_a4 );
  fGmp.a6  = fConfig->GetDoubleDef("Gmp-a6",  bba2003::kGmp_a6 );
  fGmp.a8  = fConfig->GetDoubleDef("Gmp-a8",  bba2003::kGmp_a8 );
  fGmp.a10 = fConfig->GetDoubleDef("Gmp-a10", bba2003::kGmp_a10);
  fGmp.a12 = fConfig->GetDoubleDef("Gmp-a12", bba2003::kGmp_a12);

  fGmn.a2  = fConfig->GetDoubleDef("Gmn-a2",  bba2003::kGmn_a2 );
  fGmn.a4  = fConfig->GetDoubleDef("Gmn-a4",  bba2003::kGmn_a4 );
  fGmn.a6  = fConfig->GetDoubleDef("Gmn-a6",  bba2003::kGmn_a6 );
  fGmn.a8  = fConfig->GetDoubleDef("Gmn-a8",  bba2003::kGmn_a8 );
  fGmn.a10 = fConfig->GetDoubleDef("Gmn-a10", bba2003::kGmn_a10);
  fGmn.a12 = fConfig->GetDoubleDef("Gmn-a12", bba2003::kGmn_a12);

  // Krutov parameters
  fGenA    = fConfig->GetDoubleDef("Gen-a",   bba2003::kGen_a  );
  fGenB    = fConfig->GetDoubleDef("Gen-b",   bba2003::kGen_b  );

  // Q2max
  fQ2Max   = fConfig->GetDoubleDef("Q2Max",   bba2003::kQ2Max  );
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::BBA03Fit(
                        double q2, double g0, const BBA2003Fit_t & fit) const
{
// The BBA2003 inverse polynomizal fit function for Gep,Gmp,Gmn
// Inputs:
//       q2  : momentum transfer, <0
//       g0  : G(q2=0) -> Gep=1, Gmp=mup, Gmn=mun (mu:magnetic moment)
//       fit : BBA2003 fit parameters for either Gep,Gmp,Gmn

  double Q2  = -q2;
  double Q4  =  Q2  * Q2;
  double Q6  =  Q4  * Q2;
  double Q8  =  Q6  * Q2;
  double Q10 =  Q8  * Q2;
  double Q12 =  Q10 * Q2;

  double g = g0 / (1. + fit.a2*Q2 + fit.a4*Q4 + fit.a6*Q6 +
                                   fit.a8*Q8 + fit.a10*Q10 + fit.a12*Q12);
  return g;
}
//____________________________________________________________________________

