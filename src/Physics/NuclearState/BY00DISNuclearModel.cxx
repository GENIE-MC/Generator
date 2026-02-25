 //____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 J. Tena Vidal <julia.tena-vidal@ific.uv>
 Universitat de Valencia

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/NuclearState/BY00DISNuclearModel.h"
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
  f = utils::nuclear::DISNuclFactor(x,A);

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
  // In this version, parameters are hardcoded.
  // This model is not advised, as a new version of this parameterization is available
  // in BY21DISNuclearModel
}
