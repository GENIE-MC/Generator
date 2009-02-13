//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - February 15, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Diffractive/DfrcKopeliovichPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DfrcKopeliovichPXSec::DfrcKopeliovichPXSec() :
XSecAlgorithmI("genie::DfrcKopeliovichPXSec")
{

}
//____________________________________________________________________________
DfrcKopeliovichPXSec::DfrcKopeliovichPXSec(string config) :
XSecAlgorithmI("genie::DfrcKopeliovichPXSec", config)
{

}
//____________________________________________________________________________
DfrcKopeliovichPXSec::~DfrcKopeliovichPXSec()
{

}
//____________________________________________________________________________
double DfrcKopeliovichPXSec::XSec(
        const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

//  LOG("Diffrac", pWARN)
//    << "*** No differential cross section calculation is implemented yet";

  return 1;
}
//____________________________________________________________________________
double DfrcKopeliovichPXSec::Integral(const Interaction * interaction) const
{
//  LOG("Diffrac", pWARN)
//    << "*** No diffractive cross section integration is implemented yet";

  return 1;
}
//____________________________________________________________________________
bool DfrcKopeliovichPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  if(interaction->ProcInfo().IsDiffractive()) return true;
  return false;
}
//____________________________________________________________________________
bool DfrcKopeliovichPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void DfrcKopeliovichPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DfrcKopeliovichPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DfrcKopeliovichPXSec::LoadConfig(void)
{
/*
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fGw = fConfig->GetDoubleDef("Gw", gc->GetDouble("AMNuGamma-Gw"));
*/
}
//____________________________________________________________________________

