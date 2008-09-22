//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 05, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"
#include "MEC/MECPXSec.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
MECPXSec::MECPXSec() :
XSecAlgorithmI("genie::MECPXSec")
{

}
//____________________________________________________________________________
MECPXSec::MECPXSec(string config) :
XSecAlgorithmI("genie::MECPXSec", config)
{

}
//____________________________________________________________________________
MECPXSec::~MECPXSec()
{

}
//____________________________________________________________________________
double MECPXSec::XSec(
          const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //
  // compute cross section
  //

  double xsec = 0;


  //----- The algorithm computes d^2xsec/dWdQ2
  //      Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("QMEC", pDEBUG)
     << "Jacobian for transformation to: " 
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }


  return xsec;
}
//____________________________________________________________________________
double MECPXSec::Integral(const Interaction * interaction) const
{
  //double xsec = fXSecIntegrator->Integrate(this,interaction);
  //return xsec;

  return 0;
}
//____________________________________________________________________________
bool MECPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

//  if(!proc_info.IsMEC()) return false;

  return true;
}
//____________________________________________________________________________
void MECPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MECPXSec::LoadConfig(void)
{

}
//____________________________________________________________________________

