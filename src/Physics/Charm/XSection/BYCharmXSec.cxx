#include <TMath.h>
#include <Math/Integrator.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Charm/XSection/BYCharmXSec.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
/////#include "Numerical/IntegratorI.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"

using namespace genie;
using namespace genie::constants;


//____________________________________________________________________________

double BYCharmXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const {

  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();
 // Get kinematical & init-state parameters
  const Kinematics &kinematics = interaction->Kine();
  const ProcessInfo &proc_info = interaction->ProcInfo();
  if (!proc_info.IsWeakCC()){
    return 0.;
  }

  double E = init_state.ProbeE(kRfHitNucRest);
  double ml = interaction->FSPrimLepton()->Mass();
  double Mnuc = init_state.Tgt().HitNucMass();
  double x = kinematics.x();
  double y = kinematics.y();

  double E2 = E * E;
  double ml2 = ml * ml;
  double ml4 = ml2 * ml2;
  double Mnuc2 = Mnuc * Mnuc;

  bool is_nubar_cc =
      pdg::IsAntiNeutrino(init_state.ProbePdg()) && proc_info.IsWeakCC();
  int sign = (is_nubar_cc) ? -1 : 1;

  // Calculate the DIS structure functions
  this->Calculate(interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG) << fDISSF;
#endif


  double Q2 = utils::kinematics::XYtoQ2(E, Mnuc, x, y);
  double g2 = kGF2 * kMw2 * kMw2 / TMath::Power((Q2 + kMw2), 2);
  double front_factor = (g2 * Mnuc * E) / kPi;

  // Build all dxsec/dxdy terms
  double term1 = y * (x * y + ml2 / (2 * E * Mnuc));
  double term2 = 1 - y - Mnuc * x * y / (2 * E) - ml2 / (4 * E2);
  double term3 = sign * (x * y * (1 - y / 2) - y * ml2 / (4 * Mnuc * E));
  double term4 = x * y * ml2 / (2 * Mnuc * E) + ml4 / (4 * Mnuc2 * E2);
  double term5 = -1. * ml2 / (Mnuc * E);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BYCharmXSec", pNOTICE) << "d2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2
                          << ")*F2+(" << term3 << ")*F3+(" << term4 << ")*F4+("
                          << term5 << ")*F5";
#endif
  term1 *= fF1;
  term2 *= fF2;
  term3 *= fF3;
  term4 *= fF4;
  term5 *= fF5;
  LOG("BYCharmXSec", pNOTICE) << "d2xsec/dxdy ~ "<< front_factor<<"*((" << term1 << ")*F1+(" << term2
                          << ")*F2+(" << term3 << ")*F3+(" << term4 << ")*F4+("
                          << term5 << ")*F5)";
  double xsec = front_factor * (term1 + term2 + term3 + term4 + term5);
  //LOG("BYCharmXSec", pNOTICE) << "d2xsec/dxdy = " << xsec;
  xsec = TMath::Max(xsec, 0.);

  if (kps != kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction, kPSxyfE, kps);
    xsec *= J;
  }
  return xsec;
}

void BYCharmXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYCharmXSec::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
void BYCharmXSec::LoadConfig(void)
{
  GetParam("DIS-CC-XSecScale", fCCScale);
  fDISSFModel = 0;
  fDISSFModel =
      dynamic_cast<const QPMDISStrucFuncBase *>(this->SubAlg("SFAlg"));
  assert(fDISSFModel);
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *>(this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}


//____________________________________________________________________________
bool BYCharmXSec::ValidProcess(const Interaction *interaction) const {
  if (interaction->TestBit(kISkipProcessChk))
    return true;

  const ProcessInfo &proc_info = interaction->ProcInfo();
  if (!proc_info.IsDeepInelastic())
    return false;

  const InitialState &init_state = interaction->InitState();
  int probe_pdg = init_state.ProbePdg();
  if (!pdg::IsLepton(probe_pdg))
    return false;

  if (!init_state.Tgt().HitNucIsSet())
    return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if (!pdg::IsNeutronOrProton(hitnuc_pdg))
    return false;

  return true;
}

//____________________________________________________________________________
double BYCharmXSec::Integral(const Interaction *interaction) const {
  double xsec = fXSecIntegrator->Integrate(this, interaction);
  return xsec;
}


//____________________________________________________________________________
void BYCharmXSec::Calculate(const Interaction * interaction) const
{
    // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();
  double x     = fDISSFModel->ScalingVar(interaction);
  double Q2val = fDISSFModel->Q2(interaction);
  double M = tgt.HitNucP4().M();
  double Q2pdf = TMath::Max(Q2val, fDISSFModel->fQ2min);
  // Check whether it is above charm threshold
  bool above_charm = utils::kinematics::IsAboveCharmThreshold(x, Q2val, M, fDISSFModel->fMc);
  if (!above_charm){
    return;
  }
  fDISSFModel->CalcPDFs(interaction);

  // Get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  int  nuc_pdgc    = tgt.HitNucPdg();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_p        = pdg::IsProton       ( nuc_pdgc    );
  bool is_n        = pdg::IsNeutron      ( nuc_pdgc    );
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  bool is_lepton   = pdg::IsLepton       ( probe_pdgc  );
  bool is_CC       = proc_info.IsWeakCC();

  if ( !is_nu || !is_nu) return;
  if ( !is_p && !is_n       ) return;
  if ( tgt.N() == 0 && is_n ) return;
  if ( tgt.Z() == 0 && is_p ) return;
  if ( !is_CC ) {return; }


  double switch_uv    = 0.;
  double switch_us    = 0.;
  double switch_ubar  = 0.;
  double switch_dv    = 0.;
  double switch_ds    = 0.;
  double switch_dbar  = 0.;
  double switch_s     = 0.;
  double switch_sbar  = 0.;
  double switch_c     = 0.;
  double switch_cbar  = 0.;

  if(tgt.HitQrkIsSet()) {


     int  qpdg = tgt.HitQrkPdg();
     bool sea  = tgt.HitSeaQrk();

     bool is_u    = pdg::IsUQuark     (qpdg);
     bool is_ubar = pdg::IsAntiUQuark (qpdg);
     bool is_d    = pdg::IsDQuark     (qpdg);
     bool is_dbar = pdg::IsAntiDQuark (qpdg);
     bool is_s    = pdg::IsSQuark     (qpdg);
     bool is_sbar = pdg::IsAntiSQuark (qpdg);
     bool is_c    = pdg::IsCQuark     (qpdg);
     bool is_cbar = pdg::IsAntiCQuark (qpdg);

     if      (!sea && is_u   ) { switch_uv   = 1; }
     else if ( sea && is_u   ) { switch_us   = 1; }
     else if ( sea && is_ubar) { switch_ubar = 1; }
     else if (!sea && is_d   ) { switch_dv   = 1; }
     else if ( sea && is_d   ) { switch_ds   = 1; }
     else if ( sea && is_dbar) { switch_dbar = 1; }
     else if ( sea && is_s   ) { switch_s    = 1; }
     else if ( sea && is_sbar) { switch_sbar = 1; }
     else if ( sea && is_c   ) { switch_c    = 1; }
     else if ( sea && is_cbar) { switch_cbar = 1; }
     else return;

     // make sure user inputs make sense
    if(is_nu    && is_CC && is_u   ) return;
    if(is_nu    && is_CC && is_c   ) return;
    if(is_nu    && is_CC && is_dbar) return;
    if(is_nu    && is_CC && is_sbar) return;
    if(is_nubar && is_CC && is_ubar) return;
    if(is_nubar && is_CC && is_cbar) return;
    if(is_nubar && is_CC && is_d   ) return;
    if(is_nubar && is_CC && is_s   ) return;
  } else {
    return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections

  fDISSFModel -> CalcPDFs (interaction);
  bool hasCharmContribution = (fDISSFModel->fdv_c * switch_dv   > 0) 
   || (fDISSFModel->fds_c * switch_ds   > 0) || (fDISSFModel->fs_c  * switch_s    > 0) 
   || (fDISSFModel->fc_c  * switch_cbar > 0) || (fDISSFModel->fc_c  * switch_c    > 0)
   || (fDISSFModel->fds_c * switch_dbar > 0) || (fDISSFModel->fs_c  * switch_sbar > 0);
  if (!hasCharmContribution){
    return;
  }

  // KCH = 1 if below charm threshold
  double KCH = 1.0;
  KCH = fDISSFModel->KCharm(interaction, fDISSFModel->fMc);  
  
  // Compute the K factors
  double kV_val_u = 1.;
  double kV_val_d = 1.;
  double kV_sea_u = 1.;
  double kV_sea_d = 1.;
  double kV_sea_s = 1.;
  // Axial Factors
  double kA_val_u = 1.;
  double kA_val_d = 1.;
  double kA_sea_u = 1.;
  double kA_sea_d = 1.;
  double kA_sea_s = 1.;

  if (is_n){ 
    fDISSFModel->KVectorFactors(interaction, kV_val_d, kV_val_u, kV_sea_d, kV_sea_u, kV_sea_s);
    fDISSFModel->KAxialFactors (interaction, kA_val_d, kA_val_u, kA_sea_d, kA_sea_u, kA_sea_s);
  } else {
    fDISSFModel->KVectorFactors(interaction, kV_val_u, kV_val_d, kV_sea_u, kV_sea_d, kV_sea_s);
    fDISSFModel->KAxialFactors (interaction, kA_val_u, kA_val_d, kA_sea_u, kA_sea_d, kA_sea_s);
  }


  double F2val=0, xF3val=0;
    if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {      
      // KA = KV for charm due to its mass
      q    = ( switch_dv * fDISSFModel->fdv_c * ( kV_val_d + kV_val_d ) 
             + switch_ds * fDISSFModel->fds_c * ( kV_sea_d + kV_sea_d ) ) * fDISSFModel->fVcd2 ;
      q   +=   switch_s  * fDISSFModel->fs_c  * ( kV_sea_s + kV_sea_s )   * fDISSFModel->fVcs2;

      qbar  = switch_cbar * fDISSFModel->fc_c * ( kV_sea_u + kV_sea_u ) * fDISSFModel->fVcd2;
      qbar += switch_cbar * fDISSFModel->fc_c * ( kV_sea_u + kV_sea_u ) * fDISSFModel->fVcs2;
    } else if (is_nubar) {
	    q    =   switch_c  * fDISSFModel->fc_c * ( kV_sea_u + kV_sea_u )   * fDISSFModel->fVcd2;
	    q   +=   switch_c  * fDISSFModel->fc_c * ( kV_sea_u + kV_sea_u )   * fDISSFModel->fVcs2;

	    qbar  = switch_dbar * fDISSFModel->fds_c * ( kV_sea_d + kV_sea_d ) * fDISSFModel->fVcd2;
	    qbar += switch_sbar * fDISSFModel->fs_c  * ( kV_sea_s + kV_sea_s ) * fDISSFModel->fVcs2;
    } 
    
    F2val  = (q+qbar);
    
    if (is_nu) { 
      q    = ( switch_dv * fDISSFModel->fdv_c * sqrt( kV_val_d * kV_val_d ) 
             + switch_ds * fDISSFModel->fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fDISSFModel->fVcd2;
      q   +=   switch_s  * fDISSFModel->fs_c  * sqrt( kV_sea_s * kV_sea_s )   * fDISSFModel->fVcs2;
      
      qbar  = switch_cbar * fDISSFModel->fc_c * sqrt( kV_sea_u * kV_sea_u ) * fDISSFModel->fVcd2;
      qbar += switch_cbar * fDISSFModel->fc_c * sqrt( kV_sea_u * kV_sea_u ) * fDISSFModel->fVcs2;
    }
    else if (is_nubar) {
	    q    = ( switch_c  * fDISSFModel->fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fDISSFModel->fVcd2;
	    q   += ( switch_c  * fDISSFModel->fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fDISSFModel->fVcs2;

	    qbar  = ( switch_dbar * fDISSFModel->fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fDISSFModel->fVcd2;
	    qbar += ( switch_sbar * fDISSFModel->fs_c  * sqrt( kV_sea_s * kV_sea_s ) ) * fDISSFModel->fVcs2;
    } 

    xF3val = 2*(q-qbar);
  }

  x     = fDISSFModel->ScalingVar(interaction);
  double r     = fDISSFModel->R         (interaction); // R ~ FL
  double H     = fDISSFModel->fIncludeH ? fDISSFModel->H(interaction) : 1;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pNOTICE) << "R(=FL/2xF1) = " << r;
#endif
  
  if(fDISSFModel->fUse2016Corrections) {
    const Kinematics & kinematics = interaction->Kine();
    double bjx = kinematics.x();

    double a = TMath::Power(bjx,2.) / TMath::Max(Q2val, fDISSFModel->fLowQ2CutoffF1F2);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);

    fF3 = H * KCH * xF3val/bjx;
    fF2 = F2val;
    fF1 = fF2 * KCH * 0.5 * c / bjx;
    fF5 = fF2 * 0.5 / bjx;           // Albright-Jarlskog relation
    fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)
  }
  else {
    double a = TMath::Power(x,2.) / TMath::Max(Q2val, fDISSFModel->fLowQ2CutoffF1F2);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
    fF3 = H * xF3val / x;
    fF2 = F2val;
    fF1 = fF2 * 0.5 * c / x;
    fF5 = fF2 * 0.5 / x;         // Albright-Jarlskog relation
    fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)
  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  fDISSFModel->Calculate(interaction);
  LOG("BYCharm", pNOTICE)
     << "F1-F5 = "
     << fF1 << " (" << fDISSFModel->F1() << "), " << fF2 << " (" << fDISSFModel->F2() 
     << "), " << fF3  << " (" << fDISSFModel->F3() <<"), " << fF4 <<  " (" << fDISSFModel->F4() << "), " 
     << fF5 << " (" << fDISSFModel->F5() << ")";
#endif
}