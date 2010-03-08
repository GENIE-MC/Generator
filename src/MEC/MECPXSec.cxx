//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 05, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "MEC/MECPXSec.h"
#include "Utils/KineUtils.h"

using namespace genie;

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

// We have no clue what the meson exchange current contribution is.
// This is a toy model and is not used in default event generation.

//  if(! this -> ValidProcess    (interaction) ) return 0.;
//  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();

  double Wdep  = TMath::Gaus(W, fMass, fWidth);
  double Q2dep = TMath::Power(1-Q2/fMq2d, -1.5);
  double norm  = 1.0;
  double xsec  = norm * Wdep * Q2dep;

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
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
  const InitialState & init_state = interaction -> InitState();  
  const Target & target = init_state.Tgt();

  double E  = init_state.ProbeE(kRfHitNucRest);
  int    A  = target.A();

  double norm = (E>fEc) ? fNorm*A : 0;
  return norm;
}
//____________________________________________________________________________
bool MECPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  if(!proc_info.IsMEC()) return false;

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
  fMq2d   = 1.5; // GeV
  fMass   = 1.1; // GeV
  fWidth  = 0.3; // GeV
  fEc     = 0.4; // GeV
  fNorm   = 0.2 * (1E-38 * units::cm2);;
}
//____________________________________________________________________________

