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

  double E = init_state.ProbeE(kRfHitNucRest);
  double ml = interaction->FSPrimLepton()->Mass();
  double Mnuc = init_state.Tgt().HitNucMass();
  double x = kinematics.x();
  double y = kinematics.y();

  double E2 = E * E;
  double ml2 = ml * ml;
  double ml4 = ml2 * ml2;
  double Mnuc2 = Mnuc * Mnuc;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BYCharmXSec", pDEBUG) << "Computing d2xsec/dxdy @ E = " << E
                          << ", x = " << x << ", y = " << y;
#endif


  bool is_nubar_cc =
      pdg::IsAntiNeutrino(init_state.ProbePdg()) && proc_info.IsWeakCC();
  int sign = (is_nubar_cc) ? -1 : 1;

  // Calculate the DIS structure functions
  fDISSF.Calculate(interaction);

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

//#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BYCharmXSec", pNOTICE) << "\nd2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2
                          << ")*F2+(" << term3 << ")*F3+(" << term4 << ")*F4+("
                          << term5 << ")*F5";
//#endif

  term1 *= fF1;
  term2 *= fF2;
  term3 *= fF3;
  term4 *= fF4;
  term5 *= fF5;

  double xsec = front_factor * (term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec, 0.);

//#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BYCharmXSec", pNOTICE) << "d2xsec/dxdy[FreeN] (E= " << E << ", x= " << x
                         << ", y= " << y << ") = " << xsec;
//#endif

  // The algorithm computes d^2xsec/dxdy
  // Check whether variable tranformation is needed
  if (kps != kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction, kPSxyfE, kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if (interaction->TestBit(kIAssumeFreeNucleon))
    return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must
  // have been included in the structure functions)
  const Target &target = init_state.Tgt();
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();
  xsec *= NNucl;
  xsec *= fCCScale;

  
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
      dynamic_cast<const DISStructureFuncModelI *>(this->SubAlg("SFAlg"));
  assert(fDISSFModel);
  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm
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

  // Get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();

  int  nuc_pdgc    = tgt.HitNucPdg();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_p        = pdg::IsProton       ( nuc_pdgc    );
  bool is_n        = pdg::IsNeutron      ( nuc_pdgc    );
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  bool is_lepton   = pdg::IsLepton       ( probe_pdgc  );
  bool is_dm       = pdg::IsDarkMatter   ( probe_pdgc  );
  bool is_CC       = proc_info.IsWeakCC();
  bool is_NC       = proc_info.IsWeakNC();
  bool is_EM       = proc_info.IsEM();
  bool is_dmi      = proc_info.IsDarkMatter();

  if ( !is_lepton && !is_dm ) return;
  if ( !is_p && !is_n       ) return;
  if ( tgt.N() == 0 && is_n ) return;
  if ( tgt.Z() == 0 && is_p ) return;

  // Flags switching on/off quark contributions so that this algorithm can be
  // used for both l + N -> l' + X, and l + q -> l' + q' level calculations

  double switch_uv    = 1.;
  double switch_us    = 1.;
  double switch_ubar  = 1.;
  double switch_dv    = 1.;
  double switch_ds    = 1.;
  double switch_dbar  = 1.;
  double switch_s     = 1.;
  double switch_sbar  = 1.;
  double switch_c     = 1.;
  double switch_cbar  = 1.;

  if(tgt.HitQrkIsSet()) {

     switch_uv    = 0.;
     switch_us    = 0.;
     switch_ubar  = 0.;
     switch_dv    = 0.;
     switch_ds    = 0.;
     switch_dbar  = 0.;
     switch_s     = 0.;
     switch_sbar  = 0.;
     switch_c     = 0.;
     switch_cbar  = 0.;

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
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections

  fDISSFModel -> CalcPDFs (interaction);

  // In the case of Charm CC DIS, we need a different treatment 
  bool hasCharmContribution = (fdv_c * switch_dv   > 0) || (fds_c * switch_ds   > 0) || (fs_c  * switch_s    > 0) || (fc_c  * switch_cbar > 0) || (fc_c  * switch_c    > 0) || (fds_c * switch_dbar > 0) || (fs_c  * switch_sbar > 0);


  // KCH = 1 if below charm threshold
  double KCH = 1.0;
  if (applyCharmCorrection) KCH = KCharm(interaction, fMc);  
  
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

  bool isN = pdg::IsNeutron (nuc_pdgc);
  if (isN){ 
    // if target is neutron swap u <-> d
    // like in this->calc
    fDISSFModel->KVectorFactors(interaction, kV_val_d, kV_val_u, kV_sea_d, kV_sea_u, kV_sea_s);
    fDISSFModel->KAxialFactors (interaction, kA_val_d, kA_val_u, kA_sea_d, kA_sea_u, kA_sea_s);
  } else {
    fDISSFModel->KVectorFactors(interaction, kV_val_u, kV_val_d, kV_sea_u, kV_sea_d, kV_sea_s);
    fDISSFModel->KAxialFactors (interaction, kA_val_u, kA_val_d, kA_sea_u, kA_sea_d, kA_sea_s);
  }

  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "K-Factors:";
  LOG("DISSF", pDEBUG) << "U: kV_val_u = " << kV_val_u << ", kV_sea_u = " << kV_sea_u;
  LOG("DISSF", pDEBUG) << "U: kA_val_u = " << kA_val_u << ", kA_sea_u = " << kA_sea_u;
  LOG("DISSF", pDEBUG) << "D: kV_val_d = " << kV_val_d << ", kV_sea_d = " << kV_sea_d;
  LOG("DISSF", pDEBUG) << "D: kA_val_d = " << kA_val_d << ", kA_sea_d = " << kA_sea_d;
#endif
  
  //
  // Compute structure functions for the EM, NC and CC cases
  //

  double F2val=0, xF3val=0;
  if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {      
      // KA = KV for charm due to its mass
      q  = ( switch_dv * fdv_c * ( kV_val_d + kV_val_d ) 
             + switch_ds * fds_c * ( kV_sea_d + kV_sea_d ) ) * fVcd2 ;
      q   +=   switch_s  * fs_c  * ( kV_sea_s + kV_sea_s )   * fVcs2;

      qbar = switch_cbar * fc_c * ( kV_sea_u + kV_sea_u ) * fVcd2;
      qbar += switch_cbar * fc_c * ( kV_sea_u + kV_sea_u ) * fVcs2;
    } else if (is_nubar) {
	    q   =   switch_c  * fc_c * ( kV_sea_u + kV_sea_u )   * fVcd2;
	    q   +=   switch_c  * fc_c * ( kV_sea_u + kV_sea_u )   * fVcs2;

	    qbar  = switch_dbar * fds_c * ( kV_sea_d + kV_sea_d ) * fVcd2;
	    qbar += switch_sbar * fs_c  * ( kV_sea_s + kV_sea_s ) * fVcs2;
    } else {
	    return;
    }
    
    F2val  = (q+qbar);
    
    if (is_nu) { 

      q   = ( switch_dv * fdv_c * sqrt( kV_val_d * kV_val_d ) 
             + switch_ds * fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fVcd2;
      q   +=   switch_s  * fs_c  * sqrt( kV_sea_s * kV_sea_s )   * fVcs2;
      

      qbar = switch_cbar * fc_c * sqrt( kV_sea_u * kV_sea_u ) * fVcd2;
      qbar += switch_cbar * fc_c * sqrt( kV_sea_u * kV_sea_u ) * fVcs2;
    }
    else if (is_nubar) {
	    q   = ( switch_c  * fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fVcd2;
	    q   += ( switch_c  * fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fVcs2;

	    qbar  = ( switch_dbar * fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fVcd2;
	    qbar += ( switch_sbar * fs_c  * sqrt( kV_sea_s * kV_sea_s ) ) * fVcs2;
    } else {
	    return;
    }

    xF3val = 2*(q-qbar);
  }



  double Q2val = this->Q2        (interaction);
  double x     = this->ScalingVar(interaction);
  double r     = this->R         (interaction); // R ~ FL
  double H     = fIncludeH ? this->H(interaction) : 1;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "R(=FL/2xF1) = " << r;
#endif
  
  if(fUse2016Corrections) {
    //It was confirmed by A.Bodek that the modified scaling variable
    //should just be used to compute the strucure functions F2 and xF3,
    //but that the usual Bjorken x should be used for the relations
    //between the structure functions.
    //For the same reason remove the freezing of Q2 at 0.8 for those relations,
    //although it has not been explicitly asked to A.Bodek if it should be done.

    const Kinematics & kinematics = interaction->Kine();
    double bjx = kinematics.x();

    double a = TMath::Power(bjx,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);

    fF3 = H * KCH * xF3val/bjx;
    fF2 = F2val;
    fF1 = fF2 * KCH * 0.5 * c / bjx;
    fF5 = fF2 * 0.5 / bjx;           // Albright-Jarlskog relation
    fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)
  }
  else {
    double a = TMath::Power(x,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
    //double a = TMath::Power(x,2.) / Q2val;
    //double c = (1. + 4. * kNucleonMass * a) / (1.+r);
    // KCH neq 1 if above charm threshold and cc interactions
    fF3 = H * xF3val / x;
    fF2 = F2val;
    fF1 = fF2 * 0.5 * c / x;
    fF5 = fF2 * 0.5 / x;         // Albright-Jarlskog relation
    fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG)
     << "F1-F5 = "
     << fF1 << ", " << fF2 << ", " << fF3 << ", " << fF4 << ", " << fF5;
#endif
}