 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/DeepInelastic/NuclearModel/JGDISNuclearModel.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________

JGDISNuclearModel::JGDISNuclearModel() :
  DISNuclearModelI("genie::JGDISNuclearModel"){}

//____________________________________________________________________________

JGDISNuclearModel::JGDISNuclearModel(string config) :
  DISNuclearModelI("genie::JGDISNuclearModel"){}

//____________________________________________________________________________

double JGDISNuclearModel::DISACorrection (const Interaction * interaction) const {
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

  if ( A >= 2 ) {
    double xv     = TMath::Min(0.65, TMath::Max(0.05, x));
    double xv2    = xv  * xv;
    double xv3    = xv2 * xv;
    double xv4    = xv3 * xv;
    double xv5    = xv4 * xv;
    double xv6    = xv5 * xv;
    double xv7    = xv6 * xv;
    double xv8    = xv7 * xv;
    double alpha = fAa0 + fAa1 * xv + fAa2 * xv2 + fAa3 * xv3 + fAa4 * xv4 + fAa5 * xv5 + fAa6 * xv6 + fAa7 * xv7 + fAa8 * xv8;
    double lnC   = fAc0 + fAc1 * TMath::Log(xv) + fAc2 * pow( TMath::Log(xv),2 );
    double C = TMath::Exp( lnC ) ;
    f *= C * pow( A, alpha ) ;
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "J. Gomez et. al. Nuclear factor for x of " << x << "  = " << f;
#endif

  return f;

}

//____________________________________________________________________________

void JGDISNuclearModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void JGDISNuclearModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void JGDISNuclearModel::LoadConfig(void)
{
  GetParam( "JG-NuclModel-c0", fAc0 ) ;
  GetParam( "JG-NuclModel-c1", fAc1 ) ;
  GetParam( "JG-NuclModel-c2", fAc2 ) ;

  GetParam( "JG-NuclModel-a0", fAa0 ) ;
  GetParam( "JG-NuclModel-a1", fAa1 ) ;
  GetParam( "JG-NuclModel-a2", fAa2 ) ;
  GetParam( "JG-NuclModel-a3", fAa3 ) ;
  GetParam( "JG-NuclModel-a4", fAa4 ) ;
  GetParam( "JG-NuclModel-a5", fAa5 ) ;
  GetParam( "JG-NuclModel-a6", fAa6 ) ;
  GetParam( "JG-NuclModel-a7", fAa7 ) ;
  GetParam( "JG-NuclModel-a8", fAa8 ) ;
}
