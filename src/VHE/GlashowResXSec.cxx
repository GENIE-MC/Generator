//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "VHE/GlashowResXSec.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GlashowResXSec::GlashowResXSec() :
XSecAlgorithmI("genie::GlashowResXSec")
{

}
//____________________________________________________________________________
GlashowResXSec::GlashowResXSec(string config) :
XSecAlgorithmI("genie::GlashowResXSec", config)
{

}
//____________________________________________________________________________
GlashowResXSec::~GlashowResXSec()
{

}
//____________________________________________________________________________
double GlashowResXSec::XSec(
          const Interaction * /*interaction*/, KinePhaseSpace_t /*kps*/) const
{
  return 0;
}
//____________________________________________________________________________
double GlashowResXSec::Integral(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial & final state information
  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);

  double me   = kElectronMass;
  double Mw   = kMw;
  double Gw   = PDGLibrary::Instance()->Find(kPdgWM)->Width();
  double Mw2  = TMath::Power(Mw,  2);
  double Mw4  = TMath::Power(Mw2, 2);
  double Gw2  = TMath::Power(Gw,  2);
  double s    = 2*me*E;
  double bw   = Mw4 / (TMath::Power(s-Mw2,2) + Gw2*Mw2);
  double xsec = kGF2*s*bw / (3*kPi);

  LOG("GlashowResXSec", pDEBUG) << "XSec (E = " << E << ") = " << xsec;
  return xsec;
}
//____________________________________________________________________________
bool GlashowResXSec::ValidProcess(const Interaction * interaction) const
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
bool GlashowResXSec::ValidKinematics(const Interaction * /*in*/) const
{
  return true;
}
//____________________________________________________________________________

