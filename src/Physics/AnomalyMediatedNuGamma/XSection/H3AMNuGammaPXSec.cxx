//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - February 15, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2008 - CB
   This cross section algorithm was added in 2.3.1. The differential cross
   section is not yet implemented - only the integrated cross section. 
   Energy dependence (6-th power) increases the cross section rapidly for
   energies > ~ nucleon mass which is outside the model validity range.
   Currently preventing that by frozing E-dependence at ~1 GeV.
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/AnomalyMediatedNuGamma/XSection/H3AMNuGammaPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
H3AMNuGammaPXSec::H3AMNuGammaPXSec() :
XSecAlgorithmI("genie::H3AMNuGammaPXSec")
{

}
//____________________________________________________________________________
H3AMNuGammaPXSec::H3AMNuGammaPXSec(string config) :
XSecAlgorithmI("genie::H3AMNuGammaPXSec", config)
{

}
//____________________________________________________________________________
H3AMNuGammaPXSec::~H3AMNuGammaPXSec()
{

}
//____________________________________________________________________________
double H3AMNuGammaPXSec::XSec(
        const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  LOG("AMNuGamma", pWARN)
    << "*** No differential cross section calculation is implemented yet";

  return 1;
}
//____________________________________________________________________________
double H3AMNuGammaPXSec::Integral(const Interaction * interaction) const
{
  // Compute the cross section for a free nucleon target
  const InitialState & init_state = interaction -> InitState();
  const Target &       target     = init_state.Tgt();
  double Ev    = init_state.ProbeE(kRfHitNucRest);

  double Ecutoff = kNucleonMass / 2;

  if(Ev > Ecutoff) return 0;

  double xsec0 = 2.2E-41 * units::cm2;
  double xsec  = xsec0 * TMath::Power(Ev,6.) * TMath::Power(0.1*fGw,4.);

  LOG("AMNuGamma", pNOTICE)
    << "*** xsec(vN->vNgamma) [free nuc](Ev="<< Ev << ") = "<< xsec;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("AMNuGamma", pDEBUG)
    << "*** xsec(vN->vNgamma) [free nuc](Ev="<< Ev << ") = "<< xsec;
#endif

  // If requested return the free xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Scale for the number of scattering centers at the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();
  xsec*=NNucl;

  return xsec;
}
//____________________________________________________________________________
bool H3AMNuGammaPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  if(interaction->ProcInfo().IsAMNuGamma()) return true;
  return false;
}
//____________________________________________________________________________
bool H3AMNuGammaPXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void H3AMNuGammaPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void H3AMNuGammaPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void H3AMNuGammaPXSec::LoadConfig(void)
{

  GetParam( "AMNuGamma-Gw", fGw ) ;

}
//____________________________________________________________________________

