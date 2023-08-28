//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/BWFunc.h"
#include "Physics/Resonance/XSection/MAIDRESVectorXSec.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectorXSec::MAIDRESVectorXSec() :
XSecAlgorithmI("genie::MAIDRESVectorXSec")
{

}
//____________________________________________________________________________
MAIDRESVectorXSec::MAIDRESVectorXSec(string config) :
XSecAlgorithmI("genie::MAIDRESVectorXSec", config)
{

}
//____________________________________________________________________________
MAIDRESVectorXSec::~MAIDRESVectorXSec()
{
}
//____________________________________________________________________________
double MAIDRESVectorXSec::XSec( const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics & kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  double W  = kinematics.W();
  double q2 = kinematics.q2();

  // Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("ReinSehgalRes", pDEBUG)
         << "RES/DIS Join Scheme: XSec[RES, W=" << W
         << " >= Wcut=" << fWcut << "] = 0";
#endif
       return 0;
    }
  }

  bool is_EM     = proc_info.IsEM();
  if( !is_EM ) return 0 ;

  double Mnuc = target.HitNucMass() ;
  double Mnuc2 = TMath::Power(Mnuc,2) ;
  double W2 = TMath::Power(W,2);
  double k = 0.5 * ( W2 - Mnuc2) / Mnuc ;
  double v = k - 0.5 * q2/Mnuc ;
  double v2 = TMath::Power( v, 2 ) ;
  double Q2 = v2 - q2 ; 
  double E = init_state.ProbeE(kRfHitNucRest) ;
  double Eprime = E - v ; 
  TLorentzVector ElectronMom = *init_state.GetProbeP4(kRfLab) ;
  TLorentzVector OutElectronMom = kinematics.FSLeptonP4() ; 
  double theta = OutElectronMom.Theta();
  double q3Vect2 = pow( ( ElectronMom - OutElectronMom ).P(),2) ; 

  double epsilon = 1 / ( 1 + 2 * ( q3Vect2 / Q2 ) * TMath::Power( tan( theta ) * 0.5, 2 ) ) ;      
  double Gamma = ( kAem * 0.5 / pow(kPi,2) ) * ( Eprime / E ) * ( k / Q2 ) / ( 1 - epsilon ) ; 
  double xsecT = XSecT(interaction);
  double xsecL = XSecL(interaction);
  
  double xsec = Gamma * ( xsecT + epsilon * xsecL ) ; 
  xsec = TMath::Max(0.,xsec) ; 

  xsec *= fRESScaling ; 

  return xsec ; 
}
//____________________________________________________________________________
double MAIDRESVectorXSec::XSecT( const Interaction * interaction ) const
{  

  // Get the electron, hit nucleon & current
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  int  nucpdgc   = target.HitNucPdg();
  bool is_p      = pdg::IsProton  (nucpdgc);

  RESVectFormFactorsI * vffmodel = 0;
  if (is_p) { vffmodel = fVFFEMp;}
  else      { vffmodel = fVFFEMn;}
  
  if( ! vffmodel ) return 0 ; 

  RESVectFFAmplitude vffampl = vffmodel->Compute(*interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RSHAmpl", pDEBUG)
     << "Helicity Amplitudes for RES = " << resname << " : " << vffampl;
#endif
  
  TLorentzVector ElectronMom = *init_state.GetProbeP4(kRfLab) ;  

  double Ampl2A12 = vffampl.Ampl2A12() ; 
  double Ampl2A32 = vffampl.Ampl2A32() ; 

  double factor = GetXSecFactors(interaction);
  double xsecT = factor * ( Ampl2A12 + Ampl2A32 ) ; 
  return xsecT ; 
}
//____________________________________________________________________________
double MAIDRESVectorXSec::XSecL( const Interaction * interaction ) const
{

  // Get the electron, hit nucleon & current
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  int  nucpdgc   = target.HitNucPdg();
  bool is_p      = pdg::IsProton  (nucpdgc);

  RESVectFormFactorsI * vffmodel = 0;
  if (is_p) { vffmodel = fVFFEMp;}
  else      { vffmodel = fVFFEMn;}
  
  if( ! vffmodel ) return 0 ; 

  RESVectFFAmplitude vffampl = vffmodel->Compute(*interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RSHAmpl", pDEBUG)
     << "Helicity Amplitudes for RES = " << resname << " : " << vffampl;
#endif

  double Ampl2S12 = vffampl.Ampl2S12() ; 

  double Mnuc = target.HitNucMass() ;
  double Mnuc2 = TMath::Power(Mnuc,2) ; 

  Resonance_t resonance = interaction->ExclTag().Resonance();
  double MR = utils::res::Mass(resonance) ;
  double MR2 = TMath::Power(MR,2);

  const Kinematics & kinematics = interaction -> Kine(); 
  double W  = kinematics.W();
  double W2 = TMath::Power(W,2);
  double k = 0.5 * ( W2 - Mnuc2) / Mnuc ; 
  double q2 = kinematics.q2();
  double v = k - 0.5 * q2/Mnuc ;
  double v2 = TMath::Power( v, 2 ) ;
  double Q2 = v2 - q2 ; 

  TLorentzVector ElectronMom = *init_state.GetProbeP4(kRfLab) ;
  TLorentzVector OutElectronMom = kinematics.FSLeptonP4() ; 
  double q3Vect2 = pow( ( ElectronMom - OutElectronMom ).P(),2) ;  

  double factor = GetXSecFactors(interaction);
  double xsecL = 2 * factor * Q2 * MR2 * Ampl2S12 / ( Mnuc2 * q3Vect2 ) ; 
  return xsecL;
}
//____________________________________________________________________________
double MAIDRESVectorXSec::GetXSecFactors( const Interaction * interaction ) const
{
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  const Kinematics & kinematics = interaction -> Kine(); 
  double Mnuc = target.HitNucMass() ;
  double Mnuc2 = TMath::Power(Mnuc,2) ; 
  double W  = kinematics.W();
  double W2 = TMath::Power(W,2) ;
  double k = 0.5 * ( W2 - Mnuc2) / Mnuc ;

  Resonance_t resonance = interaction->ExclTag().Resonance();
  double MR = utils::res::Mass(resonance) ;
  double MR2 = TMath::Power(MR,2);
  double kres = 0.5 * ( W2 - MR2 ) / MR ; 
  double q2 = kinematics.q2();  
  double v = k - 0.5 * q2/Mnuc ;
  double v2 = TMath::Power( v, 2 ) ;
  double Q2 = v2 - q2 ; 
  double E = init_state.ProbeE(kRfHitNucRest) ;
  double Eprime = E - v ; 
  TLorentzVector ElectronMom = *init_state.GetProbeP4(kRfLab) ;
  TLorentzVector OutElectronMom = kinematics.FSLeptonP4() ; 
  double theta = OutElectronMom.Theta();
  double q3Vect2 = pow( ( ElectronMom - OutElectronMom ).P(),2) ; 
  double epsilon = 1 / ( 1 + 2 * ( q3Vect2 / Q2 ) * TMath::Power( tan( theta ) * 0.5, 2 ) ) ;     
 
  double Gamma = ( kAem * 0.5 / pow(kPi,2) ) * ( Eprime / E ) * ( k / Q2 ) / ( 1 - epsilon ) ; 
  double delta = MR * Gamma / ( ( pow( W2 - MR2, 2) + MR*pow(Gamma,2) ) * kPi ) ;  

  double factor = 2 * kPi * Mnuc * ( kres / k ) * delta ; 
  return factor;
}
//____________________________________________________________________________
double MAIDRESVectorXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool MAIDRESVectorXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsResonant()) return false;
  if(!xcls.KnownResonance())  return false;

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));  

  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_em   = proc_info.IsEM();
  bool l_em    = (pdg::IsChargedLepton(probe) && is_em  );

  if ( !l_em) return false;

  return true ; 
}
//____________________________________________________________________________
void MAIDRESVectorXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectorXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectorXSec::LoadConfig(void)
{
  bool good_config = true ; 
  fVFFEMp = 0 ; 
  fVFFEMn = 0 ; 

  AlgFactory * algf = AlgFactory::Instance();
  fVFFEMp  = dynamic_cast<RESVectFormFactorsI*> ( algf->AdoptAlgorithm("genie::MAIDRESVectFormFactorsEMp","Default") );
  fVFFEMn  = dynamic_cast<RESVectFormFactorsI*> ( algf->AdoptAlgorithm("genie::MAIDRESVectFormFactorsEMn","Default") );
  
  if( !fVFFEMp || fVFFEMn ) { 
    good_config = false ; 
    LOG("MAIDRESVectorXSec", pERROR) << "Failed to configure Vector Form Factor" ;
  }

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;
  this->GetParam( "RES-EM-XSecScale", fRESScaling ) ; 

  fWcut = 999999;
  if(fUsingDisResJoin) {
    this->GetParam( "Wcut", fWcut ) ;
  }

  // NeuGEN limits in the allowed resonance phase space:
  // W < min{ Wmin(physical), (res mass) + x * (res width) }
  // It limits the integration area around the peak and avoids the
  // problem with huge xsec increase at low Q2 and high W.
  // In correspondence with Hugh, Rein said that the underlying problem
  // are unphysical assumptions in the model.
  this->GetParamDef( "MaxNWidthForN2Res", fN2ResMaxNWidths, 2.0 ) ;
  this->GetParamDef( "MaxNWidthForN0Res", fN0ResMaxNWidths, 6.0 ) ;
  this->GetParamDef( "MaxNWidthForGNRes", fGnResMaxNWidths, 4.0 ) ;

  this->GetParam("FermiMomentumTable", fKFTable);
  this->GetParam("RFG-UseParametrization", fUseRFGParametrization);
  this->GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

  // load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  if( !fXSecIntegrator ) {
    good_config = false ;
    LOG("MAIDRESVectorXSec", pERROR) << "XSec Integrator is not initialized" ;
  }

  if( ! good_config ) {
    LOG("MAIDRESVectorXSec", pERROR) << "Configuration has failed.";
    exit(78) ;
  }
}
//____________________________________________________________________________
