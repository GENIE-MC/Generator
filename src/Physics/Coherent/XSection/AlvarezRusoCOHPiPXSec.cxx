//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Dennis
         University of Warwick, Rutherford Appleton Laboratory,
         October 5, 2012

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/AlvarezRusoCOHPiPXSec.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

#include "Physics/Coherent/XSection/ARConstants.h"
#include "Physics/Coherent/XSection/ARSampledNucleus.h"
#include "Physics/Coherent/XSection/AlvarezRusoCOHPiPDXSec.h"
#include "Physics/Coherent/XSection/AREikonalSolution.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

using namespace alvarezruso;

//____________________________________________________________________________
AlvarezRusoCOHPiPXSec::AlvarezRusoCOHPiPXSec() :
XSecAlgorithmI("genie::AlvarezRusoCOHPiPXSec")
{
  fMultidiff = NULL;
  fLastInteraction = NULL;
}
//____________________________________________________________________________
AlvarezRusoCOHPiPXSec::AlvarezRusoCOHPiPXSec(string config) :
XSecAlgorithmI("genie::AlvarezRusoCOHPiPXSec", config)
{
  fMultidiff = NULL;
  fLastInteraction = NULL;
}
//____________________________________________________________________________
AlvarezRusoCOHPiPXSec::~AlvarezRusoCOHPiPXSec()
{
  if (fMultidiff) delete fMultidiff;
}
//____________________________________________________________________________
double AlvarezRusoCOHPiPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  
  int A       = init_state.Tgt().A(); // mass number
  int Z       = init_state.Tgt().Z(); // atomic number
  double E_nu = init_state.ProbeE(kRfLab); // neutrino energy
  
  const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  const TLorentzVector p4_pi  = kinematics.HadSystP4();
  double E_lep = p4_lep.E();
 
  if (fLastInteraction!=interaction) {
    if (fMultidiff != NULL) {
      delete fMultidiff;
      fMultidiff = NULL;
    }

    current_t current;
    if ( interaction->ProcInfo().IsWeakCC() ) {
      current = kCC;
    }
    else if ( interaction->ProcInfo().IsWeakNC() ) {
      current = kNC;
    }
    else {
      LOG("AlvarezRusoCohPi",pDEBUG)<<"Unknown current for AlvarezRuso implementation";
      return 0.;
    }
    
    flavour_t flavour;
    if ( init_state.ProbePdg() == 12 || init_state.ProbePdg() == -12) {
      flavour=kE;
    }
    else if ( init_state.ProbePdg() == 14 || init_state.ProbePdg() == -14) {
      flavour=kMu;
    }
    else if ( init_state.ProbePdg() == 16 || init_state.ProbePdg() == -16) {
      flavour=kTau;
    }
    else {
      LOG("AlvarezRusoCohPi",pDEBUG)<<"Unknown probe for AlvarezRuso implementation";
      return 0.;
    }

    nutype_t nutype;
    if ( init_state.ProbePdg() > 0) {
      nutype = kNu;
    } else {
      nutype = kAntiNu;
    }
 
    fMultidiff = new AlvarezRusoCOHPiPDXSec(Z, A ,current, flavour, nutype);
    fLastInteraction = interaction;
  }

  double xsec = fMultidiff->DXSec(E_nu, E_lep, p4_lep.Theta(), p4_lep.Phi(), p4_pi.Theta(), p4_pi.Phi());
  xsec = xsec * 1E-38 * units::cm2;
  
  if (kps != kPSElOlOpifE) {
    xsec *= utils::kinematics::Jacobian(interaction, kPSElOlOpifE, kps );
  }
  
  return (xsec);
}
//____________________________________________________________________________
double AlvarezRusoCOHPiPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlvarezRusoCOHPiPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  int nu = init_state.ProbePdg();

  if (!proc_info.IsCoherentProduction())  return false;
  if (!proc_info.IsWeak())      return false;
  if (target.HitNucIsSet())     return false;
  if (!(target.A()>1))          return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void AlvarezRusoCOHPiPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoCOHPiPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoCOHPiPXSec::LoadConfig(void)
{
  //AlgConfigPool * confp = AlgConfigPool::Instance();
  //const Registry * gc = confp->GlobalParameterList();

  /*fRo      = fConfig->GetDoubleDef("COH-Ro",            gc->GetDouble("COH-Ro"));
  fMa      = fConfig->GetDoubleDef("Ma",            gc->GetDouble("COH-Ma"));
  fa4      = fConfig->GetDoubleDef("a4",            gc->GetDouble("COHAR-a4"));
  fa5      = fConfig->GetDoubleDef("a5",            gc->GetDouble("COHAR-a5"));
  fb4      = fConfig->GetDoubleDef("b4",            gc->GetDouble("COHAR-b4"));
  fb5      = fConfig->GetDoubleDef("b5",            gc->GetDouble("COHAR-b5"));
  ffPi     = fConfig->GetDoubleDef("fPi",           gc->GetDouble("COHAR-fPi"));
  ffStar   = fConfig->GetDoubleDef("fStar",         gc->GetDouble("COHAR-fStar"));*/


  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________

