//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/BBA03ELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA03ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA03ELFormFactorsModel", config)
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
  double q2  = interaction->Kine().q2();

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
  double q2  = interaction->Kine().q2();
  double gmp = this->BBA03Fit(q2, fMuP, fGmp);
  return gmp;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double q2  = interaction->Kine().q2();

  const Target & tgt = interaction->InitState().Tgt();

  double M   = tgt.HitNucMass();             // Mnucl
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
  double q2  = interaction->Kine().q2();
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

  // BBA2003 fit coefficients
  GetParam( "BBA03-Gep-a2", fGep.a2 ) ;
  GetParam( "BBA03-Gep-a4", fGep.a4 ) ;
  GetParam( "BBA03-Gep-a6", fGep.a6 ) ;
  GetParam( "BBA03-Gep-a8", fGep.a8 ) ;
  GetParam( "BBA03-Gep-a10", fGep.a10 ) ;
  GetParam( "BBA03-Gep-a12", fGep.a12 ) ;

  GetParam( "BBA03-Gmp-a2", fGmp.a2 ) ;
  GetParam( "BBA03-Gmp-a4", fGmp.a4 ) ;
  GetParam( "BBA03-Gmp-a6", fGmp.a6 ) ;
  GetParam( "BBA03-Gmp-a8", fGmp.a8 ) ;
  GetParam( "BBA03-Gmp-a10", fGmp.a10 ) ;
  GetParam( "BBA03-Gmp-a12", fGmp.a12 ) ;

  GetParam( "BBA03-Gmn-a2", fGmn.a2 ) ;
  GetParam( "BBA03-Gmn-a4", fGmn.a4 ) ;
  GetParam( "BBA03-Gmn-a6", fGmn.a6 ) ;
  GetParam( "BBA03-Gmn-a8", fGmn.a8 ) ;
  GetParam( "BBA03-Gmn-a10", fGmn.a10 ) ;
  GetParam( "BBA03-Gmn-a12", fGmn.a12 ) ;

  // Krutov parameters
  GetParam( "BBA03-Gen-a", fGenA ) ;
  GetParam( "BBA03-Gen-b", fGenB ) ;

  // Q2max
  GetParam( "BBA03-Q2Max", fQ2Max ) ;

  // vector mass
  GetParam( "EL-Mv",fMv ) ;
  fMv2 = TMath::Power(fMv,2);

  // anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;
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

