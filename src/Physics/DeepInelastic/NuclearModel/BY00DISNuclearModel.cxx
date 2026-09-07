 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org


 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/DeepInelastic/NuclearModel/BY00DISNuclearModel.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________

BY00DISNuclearModel::BY00DISNuclearModel() :
  DISNuclearModelI("genie::BY00DISNuclearModel"){}

//____________________________________________________________________________

BY00DISNuclearModel::BY00DISNuclearModel(string config) :
  DISNuclearModelI("genie::BY00DISNuclearModel"){}

//____________________________________________________________________________

double BY00DISNuclearModel::DISACorrection (const Interaction * interaction) const {
  if ( !interaction ) return 0;
  double f = 1.;

  // Nuclear modification to Fi
  // The scaling variable can be overwritten to include corrections

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return 1.0;
  if( interaction->TestBit(kINoNuclearCorrection) ) return 1.0;

  const Target & tgt  = interaction->InitState().Tgt();

  //   The x used for computing the DIS Nuclear correction factor should be the
  //   experimental x, not the rescaled x or off-shell-rest-frame version of x
  //   (i.e. selected x).  Since we do not have access to experimental x at this
  //   point in the calculation, just use selected x.
  const Kinematics & kine  = interaction->Kine();
  double x  = kine.x();
  int    A = tgt.A();

  // Adapted from NeuGEN's nuc_factor(). Kept original comments from Hugh.

  double xv     = TMath::Min(0.75, x);
  double xv2    = xv  * xv;
  double xv3    = xv2 * xv;
  double xv4    = xv3 * xv;
  double xv5    = xv4 * xv;
  double xvp    = TMath::Power(xv, 14.417);
  double expaxv = TMath::Exp(fFef3*xv);

  // first factor goes from free nucleons to deuterium
  if(A >= 2) {
    f= f2HScale*( f2Hf0 + f2Hf1*xv + f2Hf2*xv2 + f2Hf3*xv3 + f2Hf4*xv4 + f2Hf5*xv5);
  }
  // 2nd factor goes from deuterium to iso-scalar iron
  if(A > 2) {
    f *= ( fFef0 + fFef1 * xv + fFef2 * expaxv + fFef4 * xvp ) ;
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pDEBUG) << "Nuclear factor for x of " << x << "  = " << f;
#endif

  return f;

}

//____________________________________________________________________________

void BY00DISNuclearModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BY00DISNuclearModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}

//____________________________________________________________________________

void BY00DISNuclearModel::LoadConfig(void)
{
  GetParam( "BY00-NuclModel-2H-Scale", f2HScale ) ;
  GetParam( "BY00-NuclModel-2H-f0",    f2Hf0 ) ;
  GetParam( "BY00-NuclModel-2H-f1",    f2Hf1 ) ;
  GetParam( "BY00-NuclModel-2H-f2",    f2Hf2 ) ;
  GetParam( "BY00-NuclModel-2H-f3",    f2Hf3 ) ;
  GetParam( "BY00-NuclModel-2H-f4",    f2Hf4 ) ;
  GetParam( "BY00-NuclModel-2H-f5",    f2Hf5 ) ;

  GetParam( "BY00-NuclModel-Fe-f0",    fFef0 ) ;
  GetParam( "BY00-NuclModel-Fe-f1",    fFef1 ) ;
  GetParam( "BY00-NuclModel-Fe-f2",    fFef2 ) ;
  GetParam( "BY00-NuclModel-Fe-f3",    fFef3 ) ;
  GetParam( "BY00-NuclModel-Fe-f4",    fFef4 ) ;
}
