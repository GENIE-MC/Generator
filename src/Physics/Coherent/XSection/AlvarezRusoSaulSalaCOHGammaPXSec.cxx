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

#include "Physics/Coherent/XSection/AlvarezRusoSaulSalaCOHGammaPXSec.h"
#include "Physics/Coherent/XSection/NCgamma_Diff_Cross_Section.h" 

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace NC_gamma;

//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXSec::AlvarezRusoSaulSalaCOHGammaPXSec() :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXSec")
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXSec::AlvarezRusoSaulSalaCOHGammaPXSec(string config) :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXSec", config)
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXSec::~AlvarezRusoSaulSalaCOHGammaPXSec()
{
  ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXSec::XSec( const Interaction * interaction, 
					       KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  
  int A       = init_state.Tgt().A(); // mass number
  int Z       = init_state.Tgt().Z(); // atomic number
  
  std::string mode;
  if ( init_state.ProbePdg() > 0) {
    mode = "nu";
  } else {
    mode = "nubar";
  }
  
  std::string nucleus;
  if(Z == 6 && A == 12) nucleus = "12C";
  else if(Z == 18 && A == 40) nucleus = "40Ar";
  else {
    LOG( "AlvarezRusoSaulSala",pERROR) << "WARNING: Z = " << Z 
				       <<" and A = " << A 
				       << ";  wrong nucleus for NCgamma cross section. " ;
    return -1. ; 
  }

  Diff_Cross_Section dxsec(mode, nucleus) ;

  const TLorentzVector * p4_nu = init_state.GetProbeP4( kRfLab ) ;

  const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  const TLorentzVector p4_gamma = kinematics.HadSystP4();
  
  double theta_l = p4_nu -> Angle( p4_lep.Vect() ) ;
  double theta_g = p4_nu -> Angle( p4_gamma.Vect() ) ;
  
  TVector3 vers_l = p4_nu -> Vect().Cross( p4_lep.Vect() ) ;
  TVector3 vers_g = p4_nu -> Vect().Cross( p4_gamma.Vect() ) ;

  double phi_g = vers_l.Angle( vers_g ) ;
  // this angle is also the angle between the versors of the two planes
  // the return angle is in the region [0,pi] but the differential cross section is invariant under 
  // phi_g -> - phi_g transformation so this is enough to evaluate the xsec 
  // and it saves a lot of calculations 
  
  //LOG( "AlvarezRusoSaulSala",pFATAL ) << "Phase Space. Lepton " << p4_lep << " \t Photon: " << p4_gamma ;

  double xsec = dxsec.getDiffCrossSection( p4_nu -> E(), 
					   p4_lep.E(), 
					   theta_l, 
					   theta_g, phi_g );

  // This cross section is the following
  //
  //             d4 sigma
  // ------------------------------
  //  dE_g dTheta_l dTheta_g dPhi_g
  //
  
  // The output of Diff_Cross_Section::getDiffCrossSection is in  GeV^-3 rad^-3
  // which is the desired output of GENIE as well (natural units)  

  delete p4_nu ;

  if (kps != kPSEgTlTgPgfE ) {
    xsec *= utils::kinematics::Jacobian(interaction, kPSEgTlTgPgfE, kps );
  }

  return xsec ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlvarezRusoSaulSalaCOHGammaPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();
  const XclsTag &      xcls       = interaction -> ExclTag() ;

  int nu = init_state.ProbePdg();

  if ( ! proc_info.IsCoherentProduction() )  return false;
  if ( ! proc_info.IsWeakNC() )              return false;

  if ( target.A() <= 1 )          return false;  
  if ( target.HitNucIsSet() )     return false;
  
  int A = target.A() ;
  int Z = target.Z() ;
  
  // only 12C and 40Ar are implemented at the moment 
  if ( ( A!=12 || Z!= 6 ) && ( A!=40 || Z!=18 ) ) return false ;
  
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;
  
  if ( xcls.NSingleGammas() != 1 ) return false ;

  return true;
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXSec::LoadConfig(void)
{
  
  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________

