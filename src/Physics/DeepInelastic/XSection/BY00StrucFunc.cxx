//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/DeepInelastic/XSection/BY00StrucFunc.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BY00StrucFunc::BY00StrucFunc() :
QPMDISStrucFuncBase("genie::BY00StrucFunc")
{
  this->Init();
}
//____________________________________________________________________________
BY00StrucFunc::BY00StrucFunc(string config):
QPMDISStrucFuncBase("genie::BY00StrucFunc", config)
{
  this->Init();
}
//____________________________________________________________________________
BY00StrucFunc::~BY00StrucFunc()
{

}
//____________________________________________________________________________
void BY00StrucFunc::Configure(const Registry & config)
{
// Overload Algorithm::Configure() to read the config. registry and set
// private data members.
// QPMDISStrucFuncBase::Configure() creates the owned PDF object that gets
// configured with the specified PDFModelI
// For the ReadBYParams() method see below

  QPMDISStrucFuncBase::Configure(config);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BY00StrucFunc::Configure(string param_set)
{
  QPMDISStrucFuncBase::Configure(param_set);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BY00StrucFunc::ReadBYParams(void)
{
// Get the Bodek-Yang model parameters A,B,Csea,Cv1,Cv2 from the config.
// registry and set some private data members so as not to accessing the
// registry at every calculation.
//
  GetParam( "BY-A", fA ) ;
  GetParam( "BY-B", fB ) ;
  GetParam( "BY-CsU", fCsU ) ;
  GetParam( "BY-CsD", fCsD ) ;
  GetParam( "BY-Cv1U", fCv1U ) ;
  GetParam( "BY-Cv2U", fCv2U ) ;
  GetParam( "BY-Cv1D", fCv1D ) ;
  GetParam( "BY-Cv2D", fCv2D ) ;

}
//____________________________________________________________________________
void BY00StrucFunc::Init(void)
{
  fA    = 0;
  fB    = 0;
  fCsU  = 0;
  fCsD  = 0;
  fCv1U = 0;
  fCv2U = 0;
  fCv1D = 0;
  fCv2D = 0;
}
//____________________________________________________________________________
double BY00StrucFunc::ScalingVar(const Interaction * interaction, double Mf) const
{
// Overrides QPMDISStrucFuncBase::ScalingVar() to compute the BY scaling var

  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  double myQ2 = this->Q2(interaction);
  //myQ2 = TMath::Max(Q2,fQ2min);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "Q2 at scaling var calculation = " << myQ2;
#endif
  double a  = TMath::Power( 2*kProtonMass*x, 2 ) / myQ2;
  double xw =  2*x*(myQ2+fB) / (myQ2*(1.+TMath::Sqrt(1+a)) +  2*fA*x);

  // When the final lepton is heavy, i.e. charm production, we need to use the
  // slow rescaling correction [10.1088/0954-3899/35/5/053101]
  // This is unused when Mf = 0;
  xw *= (1 + Mf * Mf / myQ2);
  return xw;
}
//____________________________________________________________________________
void BY00StrucFunc::KVectorFactors(const Interaction * interaction,
	         double & kuv, double & kdv, double & kus, double & kds, double & ks) const
{
// Overrides QPMDISStrucFuncBase::KFactors() to compute the BY K factors for
// u(valence), d(valence), u(sea), d(sea);

  double myQ2  = this->Q2(interaction);
  double GD  = 1. / TMath::Power(1.+myQ2/0.71, 2); // p elastic form factor
  double GD2 = TMath::Power(GD,2);

  kuv = (1.-GD2)*(myQ2+fCv2U)/(myQ2+fCv1U); // K - u(valence)
  kdv = (1.-GD2)*(myQ2+fCv2D)/(myQ2+fCv1D); // K - d(valence)
  kus = myQ2/(myQ2+fCsU);                   // K - u(sea)
  kds = myQ2/(myQ2+fCsD);                   // K - d(sea)
  ks =  myQ2/(myQ2+fCsD);                   // K - s(sea)
}

void BY00StrucFunc::KAxialFactors(const Interaction * interaction,
	         double & kuv, double & kdv, double & kus, double & kds, double & ks) const
{
// Overrides QPMDISStrucFuncBase::KFactors() to compute the BY K factors for
// u(valence), d(valence), u(sea), d(sea);
  KVectorFactors(interaction, kuv, kdv, kus, kds, ks);
}
//____________________________________________________________________________
