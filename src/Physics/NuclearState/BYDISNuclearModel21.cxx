 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/NuclearState/BYDISNuclearModel21.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________

BYDISNuclearModel21::BYDISNuclearModel21() :
  DISNuclearModelI("genie::BYDISNuclearModel21"){}

//____________________________________________________________________________

BYDISNuclearModel21::BYDISNuclearModel21(string config) :
  DISNuclearModelI("genie::BYDISNuclearModel21"){}

//____________________________________________________________________________

double BYDISNuclearModel21::DISACorrection (const Interaction * interaction) const {
  if ( !interaction ) return 0; 
  double f = 1.;
  
  // Nuclear modification to Fi
  // The scaling variable can be overwritten to include corrections

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return f;
  if( interaction->TestBit(kINoNuclearCorrection) ) return f;

  const Target & tgt  = interaction->InitState().Tgt();
  const Kinematics & kinematics = interaction->Kine();
  double x = kinematics.x();
  int    A = tgt.A();
  
  double xv     = TMath::Min(0.75, TMath::Max(0.05, x));
  double xv2    = xv  * xv;
  double xv3    = xv2 * xv;
  double xv4    = xv3 * xv;
  double xv5    = xv4 * xv;
  
  // first factor goes from free nucleons to deuterium
  if(A >= 2) {
    f*= 0.985*(1.+0.422*xv - 2.745*xv2 + 7.570*xv3 - 10.335*xv4 + 5.422*xv5);
  }
  
  // Computing target-mass-corrected scaling variable
  double Q2 = kinematics.Q2();
  double M  = kNucleonMass;
  double M2 = pow( M, 2);
  
  double chiTM = 2 * xv / (1+sqrt(1+4*xv2*M2/Q2) );
  // Fermi motion is already accounted for at high chiTM. To avoid double-counting, we limit the chiTM range:
  if( chiTM > 0.65 ) chiTM = 0.65;
  
  // 2nd factor goes from deuterium to iso-scalar iron
  if(A > 2) {
    f *= (1.096 - 0.38*chiTM - 0.3 * TMath::Exp(-23*chiTM) + 8 * pow(chiTM,15) ) ;								
  }
  
  if( A == 197 || A == 208 ) { // Gold nd Lead F2(Au,Pb)/F2(Fe)
    f *= (0.932 + 2.461 * chiTM - 24.23*pow(chiTM,2) + 101.03*pow(chiTM,3) - 203.47*pow(chiTM,4) + 193.85*pow(chiTM,5) - 69.82*pow(chiTM,6)); 
  }
  
  if( A == 12 ) { // F2(Fe)/F2(C)
    f /= ( 0.919 + 1.844*chiTM - 12.73*pow(chiTM,2) + 36.89*pow(chiTM,3) - 46.77*pow(chiTM,4) + 21.22*pow(chiTM,5));
  }
  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Nuclear factor for x of " << x << "  = " << f;
#endif

  return f;

}

//____________________________________________________________________________

void BYDISNuclearModel21::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BYDISNuclearModel21::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BYDISNuclearModel21::LoadConfig(void)
{
 //this->GetParam("SetName",  fSetName );
 // this->GetParam("MemberID", fMemberID);

}
