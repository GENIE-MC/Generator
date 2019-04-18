//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/GlashowResonance/XSection/GLRESPXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GLRESPXSec::GLRESPXSec() :
XSecAlgorithmI("genie::GLRESPXSec")
{

}
//____________________________________________________________________________
GLRESPXSec::GLRESPXSec(string config) :
XSecAlgorithmI("genie::GLRESPXSec", config)
{

}
//____________________________________________________________________________
GLRESPXSec::~GLRESPXSec()
{

}
//____________________________________________________________________________
double GLRESPXSec::XSec(
   const Interaction * /*interaction*/, KinePhaseSpace_t /*kps*/) const
{
  return 0;
}
//____________________________________________________________________________
double GLRESPXSec::Integral(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);

  double gf    = kGF2/(3*kPi);
  double me    = kElectronMass;
  double Mw    = kMw;
  double Gw    = PDGLibrary::Instance()->Find(kPdgWM)->Width();
  double Mw2   = TMath::Power(Mw,  2);
  double Mw4   = TMath::Power(Mw2, 2);
  double Gw2   = TMath::Power(Gw,  2);

  double s     = 2*me*E;
  double bw    = Mw4 / (TMath::Power(s-Mw2,2) + Gw2*Mw2);
  double xsec  = gf*s*bw;

  LOG("GLRES", pDEBUG) << "XSec (E = " << E << ") = " << xsec;

  return xsec;
}
//____________________________________________________________________________
bool GLRESPXSec::ValidProcess(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  bool nuok    = pdg::IsAntiNuE(init_state.ProbePdg());
  bool nucok   = !(init_state.Tgt().HitNucIsSet());
  bool ccprcok = proc_info.IsWeakCC();

  if ( !nuok    ) return false;
  if ( !nucok   ) return false;
  if ( !ccprcok ) return false;
  
  return true;
}
//____________________________________________________________________________
bool GLRESPXSec::ValidKinematics(const Interaction * /*in*/) const
{
  return true;
}
//____________________________________________________________________________

