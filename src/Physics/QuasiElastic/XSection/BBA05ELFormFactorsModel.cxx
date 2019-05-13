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

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/BBA05ELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"

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
  //-- load the BBA2005 fit coefficients
  GetParam( "BBA05-Gep-a0", fGep.a0) ;
  GetParam( "BBA05-Gep-a1", fGep.a1) ;
  GetParam( "BBA05-Gep-a2", fGep.a2) ;

  GetParam( "BBA05-Gep-b1", fGep.b1 ) ;
  GetParam( "BBA05-Gep-b2", fGep.b2 ) ;
  GetParam( "BBA05-Gep-b3", fGep.b3 ) ;
  GetParam( "BBA05-Gep-b4", fGep.b4 ) ;

  GetParam( "BBA05-Gmp-a0", fGmp.a0 ) ;
  GetParam( "BBA05-Gmp-a1", fGmp.a1 ) ;
  GetParam( "BBA05-Gmp-a2", fGmp.a2 ) ;

  GetParam( "BBA05-Gmp-b1", fGmp.b1 ) ;
  GetParam( "BBA05-Gmp-b2", fGmp.b2 ) ;
  GetParam( "BBA05-Gmp-b3", fGmp.b3 ) ;
  GetParam( "BBA05-Gmp-b4", fGmp.b4 ) ;

  GetParam( "BBA05-Gen-a0", fGen.a0 ) ;
  GetParam( "BBA05-Gen-a1", fGen.a1 ) ;
  GetParam( "BBA05-Gen-a2", fGen.a2 ) ;

  GetParam( "BBA05-Gen-b1", fGen.b1 ) ;
  GetParam( "BBA05-Gen-b2", fGen.b2 ) ;
  GetParam( "BBA05-Gen-b3", fGen.b3 ) ;
  GetParam( "BBA05-Gen-b4", fGen.b4 ) ;

  GetParam( "BBA05-Gmn-a0", fGmn.a0 ) ;
  GetParam( "BBA05-Gmn-a1", fGmn.a1 ) ;
  GetParam( "BBA05-Gmn-a2", fGmn.a2 ) ;

  GetParam( "BBA05-Gmn-b1", fGmn.b1 ) ;
  GetParam( "BBA05-Gmn-b2", fGmn.b2 ) ;
  GetParam( "BBA05-Gmn-b3", fGmn.b3 ) ;
  GetParam( "BBA05-Gmn-b4", fGmn.b4 ) ;

  //-- anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;

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
