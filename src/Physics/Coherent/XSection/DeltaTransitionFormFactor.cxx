//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (Elastic
   package)

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/Coherent/XSection/DeltaTransitionFormFactor.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

#include "Framework/Utils/StringUtils.h" 

#include "Framework/ParticleData/BaryonResUtils.h" 

using namespace genie;


//____________________________________________________________________________
DeltaTransitionFormFactor::DeltaTransitionFormFactor() :
Algorithm("genie::DeltaTransitionFormFactor")
{

}
//____________________________________________________________________________
DeltaTransitionFormFactor::DeltaTransitionFormFactor(string config) :
Algorithm("genie::DeltaTransitionFormFactor", config)
{

}
//____________________________________________________________________________
DeltaTransitionFormFactor::~DeltaTransitionFormFactor()
{

}
//____________________________________________________________________________
double DeltaTransitionFormFactor::C3V( double Q2 ) const {

  double r = sqrt( 2.0 * fKgcm0 * constants::kProtonMass * fDeltaMass / ( constants::kPi * constants::kAem * ( fMmw2+Q2) ) );

  double egcm = ( fDeltaMass2 - Q2 - pow( constants::kProtonMass, 2 ) )/ ( 2.0 * fDeltaMass ) ; 
  double qcm = sqrt(egcm*egcm + Q2);

  double Fq = 1.0 / pow(1.0 + Q2/fParam_071, 2) ; 

  double AM = fParam_03 * (1.0 + fParam_001*Q2)*exp(-fParam_023*Q2)*(qcm/ fKgcm0)*Fq;
  double a32 = constants::kSqrt3 * (  -AM ) / 2.0 ; 

  return -r * 2.0 * a32 * constants::kProtonMass * fDeltaMass / (fMpw2+Q2);

}
//____________________________________________________________________________
double DeltaTransitionFormFactor::C3VNC( double Q2 ) const {

  return C3V(Q2) * fANC ;

}
//____________________________________________________________________________
double DeltaTransitionFormFactor::C5ANC( double Q2 ) const {

  double Fd = pow( 1.0 + Q2 / fN_Delta_Ma2, -2 ) ;
  return fN_Delta_CA5_0 * Fd;

}
//____________________________________________________________________________
void DeltaTransitionFormFactor::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaTransitionFormFactor::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaTransitionFormFactor::LoadConfig(void)
{

  fDeltaMass = utils::res::Mass( kP33_1232 ) ; 
  fDeltaMass2 = pow( fDeltaMass, 2 ) ; 

  double n_Delta_Ma ;
  GetParam( "N-Delta-Ma", n_Delta_Ma ) ;
  fN_Delta_Ma2 = pow( n_Delta_Ma, 2 ) ;

  GetParam( "N-Delta-CA5_0", fN_Delta_CA5_0 ) ;

  double mN2 = pow( constants::kProtonMass, 2 ) ; 

  fKgcm0 = ( fDeltaMass2 - mN2 ) / ( 2 * fDeltaMass ) ; 
  fMpw2 = pow( fDeltaMass + constants::kProtonMass, 2 ) ;
  fMmw2 = pow( fDeltaMass - constants::kProtonMass, 2 ) ;

  double w_A ;
  GetParam( "WeinbergAngle", w_A ) ;
  
  fANC = 1. - 2 * pow( sin( w_A ), 2 ) ; 

  GetParam( "NCG-Param03" , fParam_03 ) ;
  GetParam( "NCG-Param001" , fParam_001 ) ;
  GetParam( "NCG-Param023" , fParam_023 ) ;
  GetParam( "NCG-Param071" , fParam_071 ) ;
  
  //LOG("DeltaTransitionFormFactor", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 
  
}
//____________________________________________________________________________
