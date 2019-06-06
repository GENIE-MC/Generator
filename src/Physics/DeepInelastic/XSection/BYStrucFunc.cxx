//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2009 - CA
   Renamed to BYStrucFunc from BYStructureFuncModel

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/DeepInelastic/XSection/BYStrucFunc.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStrucFunc::BYStrucFunc() :
QPMDISStrucFuncBase("genie::BYStrucFunc")
{
  this->Init();
}
//____________________________________________________________________________
BYStrucFunc::BYStrucFunc(string config):
QPMDISStrucFuncBase("genie::BYStrucFunc", config)
{
  this->Init();
}
//____________________________________________________________________________
BYStrucFunc::~BYStrucFunc()
{

}
//____________________________________________________________________________
void BYStrucFunc::Configure(const Registry & config)
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
void BYStrucFunc::Configure(string param_set)
{
  QPMDISStrucFuncBase::Configure(param_set);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStrucFunc::ReadBYParams(void)
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
void BYStrucFunc::Init(void)
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
double BYStrucFunc::ScalingVar(const Interaction * interaction) const
{
// Overrides QPMDISStrucFuncBase::ScalingVar() to compute the BY scaling var

  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  double myQ2 = this->Q2(interaction);
  //myQ2 = TMath::Max(Q2,fQ2min);
  LOG("BodekYang", pDEBUG) << "Q2 at scaling var calculation = " << myQ2;

  double a  = TMath::Power( 2*kProtonMass*x, 2 ) / myQ2;
  double xw =  2*x*(myQ2+fB) / (myQ2*(1.+TMath::Sqrt(1+a)) +  2*fA*x);
  return xw;
}
//____________________________________________________________________________
void BYStrucFunc::KFactors(const Interaction * interaction, 
	         double & kuv, double & kdv, double & kus, double & kds) const
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
}
//____________________________________________________________________________
