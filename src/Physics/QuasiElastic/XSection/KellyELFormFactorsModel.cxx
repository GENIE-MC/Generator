//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Implementation based on J.J. Kelly, Phys.Rev.C 70 (2004) 068202

 Noah Steinberg <nsteinbe \at fnal.gov>
 Fermi National Accelerator Laboratory 
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/KellyELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
KellyELFormFactorsModel::KellyELFormFactorsModel() :
ELFormFactorsModelI("genie::KellyELFormFactorsModel")
{

}
//____________________________________________________________________________
KellyELFormFactorsModel::KellyELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::KellyELFormFactorsModel", config)
{

}
//____________________________________________________________________________
KellyELFormFactorsModel::~KellyELFormFactorsModel()
{

}
//____________________________________________________________________________
double KellyELFormFactorsModel::Gep(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gep = this->KellyFit(t,fGep);
  return gep;
}
//____________________________________________________________________________
double KellyELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gmp = this->KellyFit(t,fGmp);
  gmp *= fMuP;
  return gmp;
}
//____________________________________________________________________________
double KellyELFormFactorsModel::Gen(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gen = this->KellyFit(t,fGen);
  double gvd = this->KellyVectorDipole(interaction);
  gen *= gvd;
  return gen;
}
//____________________________________________________________________________
double KellyELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  double t   = this->tau(interaction);
  double gmn = this->KellyFit(t,fGmn);
  gmn *= fMuN;
  return gmn;
}
//____________________________________________________________________________
void KellyELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KellyELFormFactorsModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void KellyELFormFactorsModel::LoadConfig(void)
{
  //-- load the Kelly fit coefficients
  GetParam( "Kelly-Gep-a0", fGep.a0) ;
  GetParam( "Kelly-Gep-a1", fGep.a1) ;
  GetParam( "Kelly-Gep-a2", fGep.a2) ;

  GetParam( "Kelly-Gep-b1", fGep.b1 ) ;
  GetParam( "Kelly-Gep-b2", fGep.b2 ) ;
  GetParam( "Kelly-Gep-b3", fGep.b3 ) ;
  GetParam( "Kelly-Gep-b4", fGep.b4 ) ;

  GetParam( "Kelly-Gmp-a0", fGmp.a0 ) ;
  GetParam( "Kelly-Gmp-a1", fGmp.a1 ) ;
  GetParam( "Kelly-Gmp-a2", fGmp.a2 ) ;

  GetParam( "Kelly-Gmp-b1", fGmp.b1 ) ;
  GetParam( "Kelly-Gmp-b2", fGmp.b2 ) ;
  GetParam( "Kelly-Gmp-b3", fGmp.b3 ) ;
  GetParam( "Kelly-Gmp-b4", fGmp.b4 ) ;

  GetParam( "Kelly-Gen-a0", fGen.a0 ) ;
  GetParam( "Kelly-Gen-a1", fGen.a1 ) ;
  GetParam( "Kelly-Gen-a2", fGen.a2 ) ;

  GetParam( "Kelly-Gen-b1", fGen.b1 ) ;
  GetParam( "Kelly-Gen-b2", fGen.b2 ) ;
  GetParam( "Kelly-Gen-b3", fGen.b3 ) ;
  GetParam( "Kelly-Gen-b4", fGen.b4 ) ;

  GetParam( "Kelly-Gmn-a0", fGmn.a0 ) ;
  GetParam( "Kelly-Gmn-a1", fGmn.a1 ) ;
  GetParam( "Kelly-Gmn-a2", fGmn.a2 ) ;

  GetParam( "Kelly-Gmn-b1", fGmn.b1 ) ;
  GetParam( "Kelly-Gmn-b2", fGmn.b2 ) ;
  GetParam( "Kelly-Gmn-b3", fGmn.b3 ) ;
  GetParam( "Kelly-Gmn-b4", fGmn.b4 ) ;

  //-- anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;

  //-- QE Vector mass
  GetParam( "QEL-Mv", fMv) ;
  fMv2 = TMath::Power(fMv,2);

}
//____________________________________________________________________________
double KellyELFormFactorsModel::KellyFit(
                                      double t, const KellyFit_t & fp) const
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
double KellyELFormFactorsModel::KellyVectorDipole(const Interaction * interaction) const
{

  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  double Gvd = 1. / TMath::Power( 1. - q2/fMv2 , 2);
  return Gvd;
}
//____________________________________________________________________________
double KellyELFormFactorsModel::tau(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  const Target &     target = interaction->InitState().Tgt();

  double q2 = kine.q2(); // momentum transfer, <0
  double M2 = TMath::Power(target.HitNucMass(),2); // Mnucl^2

  double t = -q2/(4*M2);
  return t;
}
//____________________________________________________________________________
