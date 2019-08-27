//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda
         University of Liverpool 
         August 23, 2019

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
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

#include "Physics/Coherent/XSection/AlvarezRusoSaulSalaCOHGammaPXsec.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::AlvarezRusoSaulSalaCOHGammaPXsec() :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXsec")
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::AlvarezRusoSaulSalaCOHGammaPXsec(string config) :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXsec", config)
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::~AlvarezRusoSaulSalaCOHGammaPXsec()
{
  ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXsec::XSec( const Interaction * interaction, 
					       KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  
  int A       = init_state.Tgt().A(); // mass number
  int Z       = init_state.Tgt().Z(); // atomic number
  double E_nu = init_state.ProbeE(kRfLab); // neutrino energy
  
  // const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  // const TLorentzVector p4_pi  = kinematics.HadSystP4();
  // double E_lep = p4_lep.E();
 
  // if (fLastInteraction!=interaction) {
  //   if (fMultidiff != NULL) {
  //     delete fMultidiff;
  //     fMultidiff = NULL;
  //   }

  //   current_t current;
  //   if ( interaction->ProcInfo().IsWeakCC() ) {
  //     current = kCC;
  //   }
  //   else if ( interaction->ProcInfo().IsWeakNC() ) {
  //     current = kNC;
  //   }
  //   else {
  //     LOG("AlvarezRusoCohPi",pDEBUG)<<"Unknown current for AlvarezRuso implementation";
  //     return 0.;
  //   }
    
  //   flavour_t flavour;
  //   if ( init_state.ProbePdg() == 12 || init_state.ProbePdg() == -12) {
  //     flavour=kE;
  //   }
  //   else if ( init_state.ProbePdg() == 14 || init_state.ProbePdg() == -14) {
  //     flavour=kMu;
  //   }
  //   else if ( init_state.ProbePdg() == 16 || init_state.ProbePdg() == -16) {
  //     flavour=kTau;
  //   }
  //   else {
  //     LOG("AlvarezRusoCohPi",pDEBUG)<<"Unknown probe for AlvarezRuso implementation";
  //     return 0.;
  //   }

  //   nutype_t nutype;
  //   if ( init_state.ProbePdg() > 0) {
  //     nutype = kNu;
  //   } else {
  //     nutype = kAntiNu;
  //   }
 
  //   fMultidiff = new AlvarezRusoCOHPiPDXSec(Z, A ,current, flavour, nutype);
  //   fLastInteraction = interaction;
  // }

  // double xsec = fMultidiff->DXSec(E_nu, E_lep, p4_lep.Theta(), p4_lep.Phi(), p4_pi.Theta(), p4_pi.Phi());
  // xsec = xsec * 1E-38 * units::cm2;
  
  // if (kps != kPSElOlOpifE) {
  //   xsec *= utils::kinematics::Jacobian(interaction, kPSElOlOpifE, kps );
  // }
  
  // return (xsec);
  return 0. ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXsec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlvarezRusoSaulSalaCOHGammaPXsec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();
  const XclsTag &      xcls       = interaction -> ExclTag() ;

  int nu = init_state.ProbePdg();

  if ( ! proc_info.IsCoherent() )  return false;
  if ( ! proc_info.IsWeakNC() )    return false;

  if ( target.A() <= 1 )          return false;  
  if ( target.HitNucIsSet() )     return false;
  
  int A = target.A() ;
  int Z = target.Z() ;
  
  // only 12C and 40Ar are implemented at the moment 
  if ( ( A!=12 || Z!= 6 ) && ( A!=40 || Z!=18 ) ) return false ;
  
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;
  
  return true;
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::LoadConfig(void)
{
  
  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________

