 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/NuclearState/JGDISNuclearModel.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________

JGDISNuclearModel::JGDISNuclearModel() :
  DISNuclearModelI("genie::JGDISNuclearModel"){}

//____________________________________________________________________________

JGDISNuclearModel::JGDISNuclearModel(string config) :
  DISNuclearModelI("genie::JGDISNuclearModel"){}

//____________________________________________________________________________

double JGDISNuclearModel::DISACorrection (const Interaction & interaction) const {
  double f = 1.;
  
  // Nuclear modification to Fi
  // The scaling variable can be overwritten to include corrections

  if( interaction.TestBit(kIAssumeFreeNucleon)   ) return f;
  if( interaction.TestBit(kINoNuclearCorrection) ) return f;

  const Target & tgt  = interaction.InitState().Tgt();
  const Kinematics & kinematics = interaction.Kine();
  double x = kinematics.x();
  int    A = tgt.A();
  
  double xv     = TMath::Min(0.75, TMath::Max(0.05, x));
  double xv2    = xv  * xv;
  double xv3    = xv2 * xv;
  double xv4    = xv3 * xv;
  double xv5    = xv4 * xv;
  double xv6    = xv5 * xv;
  double xv7    = xv6 * xv;
  double xv8    = xv7 * xv;
  double alpha = -0.070 + 2.189 * xv - 24.667 * xv2 + 145.291 * xv3 - 497.237 * xv4 + 1013.129 * xv5 - 1208.393 * xv6 +775.767 * xv7 - 205.872 * xv8; 
  double lnC   = 0.017 +0.018 * TMath::Log(xv) + 0.005 * pow( TMath::Log(xv),2 );
  double C = TMath::Exp( lnC ) ;
  f *= C * pow( A, alpha ) ;
    
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
 //this->GetParam("SetName",  fSetName );
 // this->GetParam("MemberID", fMemberID);

}
