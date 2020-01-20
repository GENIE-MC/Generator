//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/Coherent/XSection/DeltaInMediumCorrections.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

#include "Framework/Utils/StringUtils.h" 

using namespace genie;


//____________________________________________________________________________
DeltaInMediumCorrections::DeltaInMediumCorrections() :
Algorithm("genie::DeltaInMediumCorrections")
{

}
//____________________________________________________________________________
DeltaInMediumCorrections::DeltaInMediumCorrections(string config) :
Algorithm("genie::DeltaInMediumCorrections", config)
{

}
//____________________________________________________________________________
DeltaInMediumCorrections::~DeltaInMediumCorrections()
{

}
//____________________________________________________________________________
void DeltaInMediumCorrections::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaInMediumCorrections::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaInMediumCorrections::LoadConfig(void)
{
  
  //  LOG("DeltaInMediumCorrections", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 
  
}
//____________________________________________________________________________
