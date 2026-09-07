 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/DeepInelastic/NuclearModel/BY21DISNuclearModel.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________

BY21DISNuclearModel::BY21DISNuclearModel() :
  DISNuclearModelI("genie::BY21DISNuclearModel"){}

//____________________________________________________________________________

BY21DISNuclearModel::BY21DISNuclearModel(string config) :
  DISNuclearModelI("genie::BY21DISNuclearModel"){}

//____________________________________________________________________________

double BY21DISNuclearModel::DISACorrection (const Interaction * interaction) const {
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
    f*= f2HScale*( f2Hf0 + f2Hf1*xv + f2Hf2*xv2 + f2Hf3*xv3 + f2Hf4*xv4 + f2Hf5*xv5);
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
    f *= ( fFef0 + fFef1 * chiTM + fFef2 * TMath::Exp( fFef3*chiTM ) + fFef4 * pow(chiTM,15) ) ;
  }

  if( A == 197 || A == 208 ) { // Gold nd Lead F2(Au,Pb)/F2(Fe)
    f *= ( fAuf0 + fAuf1 * chiTM + fAuf2 *pow(chiTM,2) + fAuf3*pow(chiTM,3) + fAuf4*pow(chiTM,4) + fAuf5*pow(chiTM,5) + fAuf6*pow(chiTM,6));
  }

  if( A == 12 ) { // F2(Fe)/F2(C)
    f /= ( fCf0 + fCf1*chiTM + fCf2*pow(chiTM,2) + fCf3*pow(chiTM,3) + fCf4*pow(chiTM,4) + fCf5*pow(chiTM,5));
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Nuclear factor for x of " << x << "  = " << f;
#endif

  return f;

}

//____________________________________________________________________________

void BY21DISNuclearModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BY21DISNuclearModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BY21DISNuclearModel::LoadConfig(void)
{

  GetParam( "BY21-NuclModel-2H-Scale", f2HScale ) ;
  GetParam( "BY21-NuclModel-2H-f0",    f2Hf0 ) ;
  GetParam( "BY21-NuclModel-2H-f1",    f2Hf1 ) ;
  GetParam( "BY21-NuclModel-2H-f2",    f2Hf2 ) ;
  GetParam( "BY21-NuclModel-2H-f3",    f2Hf3 ) ;
  GetParam( "BY21-NuclModel-2H-f4",    f2Hf4 ) ;
  GetParam( "BY21-NuclModel-2H-f5",    f2Hf5 ) ;

  GetParam( "BY21-NuclModel-Fe-f0",    fFef0 ) ;
  GetParam( "BY21-NuclModel-Fe-f1",    fFef1 ) ;
  GetParam( "BY21-NuclModel-Fe-f2",    fFef2 ) ;
  GetParam( "BY21-NuclModel-Fe-f3",    fFef3 ) ;
  GetParam( "BY21-NuclModel-Fe-f4",    fFef4 ) ;

  GetParam( "BY21-NuclModel-C-f0",    fCf0 ) ;
  GetParam( "BY21-NuclModel-C-f1",    fCf1 ) ;
  GetParam( "BY21-NuclModel-C-f2",    fCf2 ) ;
  GetParam( "BY21-NuclModel-C-f3",    fCf3 ) ;
  GetParam( "BY21-NuclModel-C-f4",    fCf4 ) ;
  GetParam( "BY21-NuclModel-C-f5",    fCf5 ) ;

  GetParam( "BY21-NuclModel-Au-f0",    fAuf0 ) ;
  GetParam( "BY21-NuclModel-Au-f1",    fAuf1 ) ;
  GetParam( "BY21-NuclModel-Au-f2",    fAuf2 ) ;
  GetParam( "BY21-NuclModel-Au-f3",    fAuf3 ) ;
  GetParam( "BY21-NuclModel-Au-f4",    fAuf4 ) ;
  GetParam( "BY21-NuclModel-Au-f5",    fAuf5 ) ;
  GetParam( "BY21-NuclModel-Au-f6",    fAuf6 ) ;

}
