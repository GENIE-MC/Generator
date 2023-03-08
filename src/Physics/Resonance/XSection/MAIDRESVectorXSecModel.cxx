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
#include "Physics/Resonance/XSection/MAIDRESVectorXSecModel.h"
#include "Physics/Resonance/XSection/RSHelicityAmplModelI.h"
#include "Physics/Resonance/XSection/RSHelicityAmpl.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
MAIDRESVectorXSecModel::MAIDRESVectorXSecModel() :
XSecAlgorithmI("genie::MAIDRESVectorXSecModel")
{

}
//____________________________________________________________________________
MAIDRESVectorXSecModel::MAIDRESVectorXSecModel(string config) :
XSecAlgorithmI("genie::MAIDRESVectorXSecModel", config)
{

}
//____________________________________________________________________________
MAIDRESVectorXSecModel::~MAIDRESVectorXSecModel()
{
}
//____________________________________________________________________________
double MAIDRESVectorXSecModel::XSec( const Interaction * interaction, KinePhaseSpace_t kps) const
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

  // Get the input baryon resonance
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  
  // Get the electron, hit nucleon & current
  int  nucpdgc   = target.HitNucPdg();
  bool is_p      = pdg::IsProton  (nucpdgc);
  bool is_EM     = proc_info.IsEM();

  if( !is_EM ) return 0 ;

  const RESVectFormFactorsI * vffmodel = 0;
  if(is_EM) {
    if (is_p) { vffmodel = fVFFModelEMp;}
    else      { vffmodel = fVFFModelEMn;}
  }
  if( ! vffmodel ) return 0 ; 

  RESVectFFAmplitude vffampl = vffmodel->Compute(*interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RSHAmpl", pDEBUG)
     << "Helicity Amplitudes for RES = " << resname << " : " << vffampl;
#endif
  
  TLorentzVector * tempElectron = init_state.GetProbeP4(kRfLab) ; 
  TLorentzVector ElectronMom = * tempElectron ; 
  delete tempElectron ; 
  TLorentzVector OutElectronMom = kinematics.FSLeptonP4() ; 

  double Mnuc = target.HitNucMass() ;
  double Mnuc2 = TMath::Power(Mnuc,2) ; 
  double W2 = TMath::Power(W,2) ;
  double k = 0.5 * ( W2 - Mnuc2) / Mnuc ;
  double MR = utils::res::Mass(resonance) ;
  double MR2 = TMath::Power(MR,2);
  double kres = 0.5 * ( W2 - MR2 ) / MR ; 
  double v = k - 0.5 * q2/Mnuc ;
  double v2 = TMath::Power( v, 2 ) ;
  double Q2 = v2 - q2 ; 
  double E = init_state.ProbeE(kRfHitNucRest) ;
  double Eprime = E - v ; 
  double theta = OutElectronMom.Theta();
  double q3Vect2 = pow( ( ElectronMom - OutElectronMom ).P(),2) ; 
  double epsilon = 1 / ( 1 + 2 * ( q3Vect2 / Q2 ) * TMath::Power( tan( theta ) * 0.5, 2 ) ) ;     
 
  double Gamma = ( kAem * 0.5 / pow(kPi,2) ) * ( Eprime / E ) * ( k / Q2 ) / ( 1 - epsilon ) ; 
  double delta = MR * Gamma / ( ( pow( W2 - MR2, 2) + MR*pow(Gamma,2) ) * kPi ) ;  
  double Ampl2A12 = vffampl.Ampl2A12() ; 
  double Ampl2A32 = vffampl.Ampl2A32() ; 
  double Ampl2S12 = vffampl.Ampl2S12() ; 

  double xsecT = 2 * Mnuc * kPi * ( kres / k ) * ( Ampl2A12 + Ampl2A32 ) * delta ; 
  double xsecL = 4 * kPi * Mnuc * kres * Q2 * MR2 * Ampl2S12 * delta / ( Mnuc2 * k * q3Vect2 ) ; 
  
  double xsec = Gamma * ( xsecT + epsilon * xsecL ) ; 
  xsec = TMath::Max(0.,xsec) ; 

  // The algorithm computes d2xsec/dEdQ2
  // Check whether variable transformation is needed
  if( kps!=kPSQ2fE ) { 
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED_
    LOG("MAIDRESVectorXSecModel", pDEBUG) << " Jacobian transformation to: " 
					  << KinePhaseSpace::AsString(kps) << ", J = " << J ; 
#endif 
    xsec *= J ; 
  }

  // If requested, return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  int Z = target.Z() ; 
  int N = target.N() ; 
  int A = target.A() ; 
  // Take into account the number of nucleons/tgt
  int NNucl = (pdg::IsProton(nucpdgc)) ? Z : N ;
  xsec *= NNucl;

  if (fUsePauliBlocking && A != 1)
    {
      // Calculation of Pauli blocking according references:
      //
      //     [1] S.L. Adler,  S. Nussinov,  and  E.A.  Paschos,  "Nuclear
      //         charge exchange corrections to leptonic pion  production
      //         in  the (3,3) resonance  region,"  Phys. Rev. D 9 (1974)
      //         2125-2143 [Erratum Phys. Rev. D 10 (1974) 1669].
      //     [2] J.Y. Yu, "Neutrino interactions and  nuclear  effects in
      //         oscillation experiments and the  nonperturbative disper-
      //         sive  sector in strong (quasi-)abelian  fields,"  Ph. D.
      //         Thesis, Dortmund U., Dortmund, 2002 (unpublished).
      //     [3] E.A. Paschos, J.Y. Yu,  and  M. Sakuda,  "Neutrino  pro-
      //         duction  of  resonances,"  Phys. Rev. D 69 (2004) 014013
      //         [arXiv: hep-ph/0308130].

      double P_Fermi = 0.0;

      // Maximum value of Fermi momentum of target nucleon (GeV)
      if (A<6 || !fUseRFGParametrization)
	{
	  // Look up the Fermi momentum for this target
	  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
	  const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
	  P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
	}
      else {
        // Define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // Correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
      }

      double FactorPauli_RES = 1.0;

      double k0 = 0., q = 0., q0 = 0.;

      if (P_Fermi > 0.)
	{
	  k0 = (W2-Mnuc2-Q2)/(2*W);
	  k = TMath::Sqrt(k0*k0+Q2);                  // previous value of k is overridden
	  q0 = (W2-Mnuc2+kPionMass2)/(2*W);
	  q = TMath::Sqrt(q0*q0-kPionMass2);
	}

      if (2*P_Fermi < k-q)
        FactorPauli_RES = 1.0;
      if (2*P_Fermi >= k+q)
        FactorPauli_RES = ((3*k*k+q*q)/(2*P_Fermi)-(5*TMath::Power(k,4)+TMath::Power(q,4)+10*k*k*q*q)/(40*TMath::Power(P_Fermi,3)))/(2*k);
      if (2*P_Fermi >= k-q && 2*P_Fermi <= k+q)
        FactorPauli_RES = ((q+k)*(q+k)-4*P_Fermi*P_Fermi/5-TMath::Power(k-q, 3)/(2*P_Fermi)+TMath::Power(k-q, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*q*k);

      xsec *= FactorPauli_RES;
    }

  // Include scaing

  return xsec ; 

}
//____________________________________________________________________________
double MAIDRESVectorXSecModel::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool MAIDRESVectorXSecModel::ValidProcess(const Interaction * interaction) const
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
void MAIDRESVectorXSecModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectorXSecModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MAIDRESVectorXSecModel::LoadConfig(void)
{
  bool good_config = true ; 
  fVFFModelEMp = 0 ; 
  fVFFModelEMn = 0 ; 

  AlgFactory * algf = AlgFactory::Instance();
  fVFFModelEMp  = dynamic_cast<const RESVectFormFactorsI*> ( algf->GetAlgorithm("genie::MAIDRESVectFormFactorsEMp","Default") );
  fVFFModelEMn  = dynamic_cast<const RESVectFormFactorsI*> ( algf->GetAlgorithm("genie::MAIDRESVectFormFactorsEMn","Default") );
  
  if( !fVFFModelEMp || fVFFModelEMn ) { 
    good_config = false ; 
    LOG("MAIDRESVectorXSecModel", pERROR) << "Failed to configure Vector Form Factor" ;
  }

  // Use algorithm within a DIS/RES join scheme. If yes get Wcut
  this->GetParam( "UseDRJoinScheme", fUsingDisResJoin ) ;
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
    LOG("MAIDRESVectorXSecModel", pERROR) << "XSec Integrator is not initialized" ;
  }

  if( ! good_config ) {
    LOG("MAIDRESVectorXSecModel", pERROR) << "Configuration has failed.";
    exit(78) ;
  }
}
//____________________________________________________________________________
