//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   First included in v2.5.1.

*/
//____________________________________________________________________________

//#include <TMath.h>

//#include "Algorithm/AlgConfigPool.h"
//#include "Base/XSecIntegratorI.h"
//#include "Conventions/Constants.h"
//#include "Conventions/RefFrame.h"
//#include "Conventions/KineVar.h"
#include "Elastic/RosenbluthPXSec.h"
//#include "Messenger/Messenger.h"
//#include "PDG/PDGUtils.h"
//#include "Utils/KineUtils.h"
//#include "Utils/MathUtils.h"
//#include "Utils/NuclearUtils.h"

using namespace genie;
//using namespace genie::utils;
//using namespace genie::constants;

//____________________________________________________________________________
RosenbluthPXSec::RosenbluthPXSec() :
XSecAlgorithmI("genie::RosenbluthPXSec")
{

}
//____________________________________________________________________________
RosenbluthPXSec::RosenbluthPXSec(string config) :
XSecAlgorithmI("genie::RosenbluthPXSec", config)
{

}
//____________________________________________________________________________
RosenbluthPXSec::~RosenbluthPXSec()
{

}
//____________________________________________________________________________
double RosenbluthPXSec::XSec(
    const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  return 1;
}
//____________________________________________________________________________
double RosenbluthPXSec::Integral(const Interaction * /*interaction*/) const
{
  return 1.;
}
//____________________________________________________________________________
bool RosenbluthPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
void RosenbluthPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RosenbluthPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RosenbluthPXSec::LoadConfig(void)
{

}
//____________________________________________________________________________
