//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 This GENIE code was adapted from the neugen3 code co-authored by Donna Naples
 (Pittsburgh U.), Hugh Gallagher (Tufts U), and Costas Andreopoulos (RAL)

 A fix was installed (Aug 12, 2014) by Brian Tice (Rochester) so that
 the nuclear modification to the pdf should be calculated in terms
 of the experimental x, not the rescaled x. The same goes for R(x,Q2).

 A fix of the scaling variable used for the relations between structure
 functions was installed by C. Bronner and J. Morrison Jun 06, 2016
 after it was confirmed by A. Bodek that x and not the modified
 scaling variable should be used there.

 Changes required to implement the GENIE Boosted Dark Matter module
 were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PhysUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QPMDISStrucFuncBase::QPMDISStrucFuncBase() :
DISStructureFuncModelI()
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDISStrucFuncBase::QPMDISStrucFuncBase(string name) :
DISStructureFuncModelI(name)
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDISStrucFuncBase::QPMDISStrucFuncBase(string name, string config):
DISStructureFuncModelI(name, config)
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDISStrucFuncBase::~QPMDISStrucFuncBase()
{
  delete fPDF;
  delete fPDFc;
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::LoadConfig(void)
{
  LOG("DISSF", pDEBUG) << "Loading configuration...";

  //-- pdf
  const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (this->SubAlg("PDF-Set"));
  fPDF  -> SetModel(pdf_model);
  fPDFc -> SetModel(pdf_model);

  //-- get CKM elements
  GetParam( "CKM-Vcd", fVcd ) ;
  GetParam( "CKM-Vcs", fVcs ) ;
  GetParam( "CKM-Vud", fVud ) ;
  GetParam( "CKM-Vus", fVus ) ;

  fVcd2 = TMath::Power( fVcd, 2 );
  fVcs2 = TMath::Power( fVcs, 2 );
  fVud2 = TMath::Power( fVud, 2 );
  fVus2 = TMath::Power( fVus, 2 );

  //-- charm mass
  GetParam( "Charm-Mass", fMc ) ;

  //-- min Q2 for PDF evaluation
  GetParam( "PDF-Q2min", fQ2min ) ;

  //-- include R (~FL)?
  GetParam( "IncludeR", fIncludeR ) ;

  //-- include nuclear factor (shadowing / anti-shadowing / ...)
  GetParam( "IncludeNuclMod", fIncludeNuclMod ) ;

  //-- Use 2016 SF relation corrections
  GetParam( "Use2016Corrections", fUse2016Corrections ) ;

  //-- Set min for relation between 2xF1 and F2
  GetParam( "LowQ2CutoffF1F2", fLowQ2CutoffF1F2 ) ;

  //-- turn charm production off?
  GetParamDef( "Charm-Prod-Off", fCharmOff, false ) ;

  //-- include H?
  GetParam( "IncludeH", fIncludeH, false ) ;

  //-- weinberg angle
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin2thw = TMath::Power(TMath::Sin(thw), 2);


  fDISNuclCorr = dynamic_cast<const DISNuclearModelI*>(this->SubAlg("DISNuclModel"));

  LOG("DISSF", pDEBUG) << "Done loading configuration";
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::InitPDF(void)
{
                     // evaluated at:
  fPDF  = new PDF(); //   x = computed (+/-corrections) scaling var, Q2
  fPDFc = new PDF(); //   x = computed charm slow re-scaling var,    Q2
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::Calculate(const Interaction * interaction) const
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
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections
  this -> CalcPDFs (interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  fPDFc->Print(std::cout);
  LOG("DISSF", pNOTICE) << "Contribution: " 
   << " dv "   << fdv_c << " " << switch_dv   << " "
   << " ds "   << fds_c << " " << switch_ds   << " " 
   << " s "    << fs_c  << " " << switch_s    << " " 
   << " cbar " << fc_c  << " " << switch_cbar << " " 
   << " c "    << fc_c  << " " << switch_c    << " "
   << " dbar " << fds_c << " " << switch_dbar << " " 
   << " sbar " << fs_c  << " " << switch_sbar << " ";

#endif
  // In the case of Charm CC DIS, we need a different treatment 
  bool hasCharmContribution = (fdv_c * switch_dv   > 0) 
   || (fds_c * switch_ds   > 0) || (fs_c  * switch_s    > 0) 
   || (fc_c  * switch_cbar > 0) || (fc_c  * switch_c    > 0)
   || (fds_c * switch_dbar > 0) || (fs_c  * switch_sbar > 0);
  const bool applyCharmCorrection = is_CC && hasCharmContribution && !fCharmOff ;

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
    this->KVectorFactors(interaction, kV_val_d, kV_val_u, kV_sea_d, kV_sea_u, kV_sea_s);
    this->KAxialFactors (interaction, kA_val_d, kA_val_u, kA_sea_d, kA_sea_u, kA_sea_s);
  } else {
    this->KVectorFactors(interaction, kV_val_u, kV_val_d, kV_sea_u, kV_sea_d, kV_sea_s);
    this->KAxialFactors (interaction, kA_val_u, kA_val_d, kA_sea_u, kA_sea_d, kA_sea_s);
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
  
  // ***  NEUTRAL CURRENT

  // Include DM in NC
  if(is_NC || is_dmi) {

    if(!is_nu && !is_nubar && !is_dm) return;

    double GL   = (is_nu) ? ( 0.5 - (2./3.)*fSin2thw) : (     - (2./3.)*fSin2thw); // clu
    double GR   = (is_nu) ? (     - (2./3.)*fSin2thw) : ( 0.5 - (2./3.)*fSin2thw); // cru
    double GLp  = (is_nu) ? (-0.5 + (1./3.)*fSin2thw) : (       (1./3.)*fSin2thw); // cld
    double GRp  = (is_nu) ? (       (1./3.)*fSin2thw) : (-0.5 + (1./3.)*fSin2thw); // crd
    // Set the couplings to up and down quarks to be axial for DM
    if (is_dm) {
      GL  = -1.;
      GR  =  1.;
      GLp = -1.;
      GRp =  1.;
    }
    double gvu  = GL  + GR;
    double gau  = GL  - GR;
    double gvd  = GLp + GRp;
    double gad  = GLp - GRp;
    double gvu2 = TMath::Power(gvu, 2.);
    double gau2 = TMath::Power(gau, 2.);
    double gvd2 = TMath::Power(gvd, 2.);
    double gad2 = TMath::Power(gad, 2.);

    double q2   = switch_uv   * fuv * (kV_val_u * gvu2 + kA_val_u * gau2);
    q2         += (switch_us  * fus + switch_c * fc )  * (kV_sea_u * gvu2 + kA_sea_u * gau2);
    q2         += switch_dv   * fdv * ( kV_val_d * gvd2+ kA_val_d * gad2);
    q2         += switch_ds   * fds * ( kV_sea_d * gvd2+ kA_sea_d * gad2);
    q2         += switch_s    * fs  * ( kV_sea_s * gvd2+ kA_sea_s * gad2);

    double qb2  = switch_ubar * fus * ( kV_sea_u * gvu2 + kA_sea_u * gau2);
    qb2        += switch_cbar * fc  * ( kV_sea_u * gvu2 + kA_sea_u * gau2);
    qb2        += switch_dbar * fds * ( kV_sea_d * gvd2 + kA_sea_d * gad2);
    qb2        += switch_sbar * fs  * ( kV_sea_s * gvd2 + kA_sea_s * gad2);
    
    double q3   = switch_uv * sqrt( kV_val_u * kA_val_u ) * fuv * (gvu*gau);
    q3         += switch_us * sqrt( kV_sea_u * kA_sea_u ) * fus * (gvu*gau);
    q3         += switch_c  * sqrt( kV_sea_u * kA_sea_u ) * fc  * (gvu*gau);
    q3         += switch_dv * sqrt( kV_val_d * kA_val_d ) * fdv * (gvd*gad);
    q3         += switch_ds * sqrt( kV_sea_d * kA_sea_d ) * fds * (gvd*gad);
    q3         += switch_s  * sqrt( kV_sea_s * kA_sea_s ) * fs  * (gvd*gad);
   

    double qb3  = switch_ubar * sqrt( kV_sea_u * kA_sea_u )  * fus * (gvu*gau);
    qb3        += switch_cbar * sqrt( kV_sea_u * kA_sea_u )  * fc  * (gvu*gau);
    qb3        += switch_dbar * sqrt( kV_sea_d * kA_sea_d )  * fds * (gvd*gad);
    qb3        += switch_sbar * sqrt( kV_sea_s * kA_sea_s )  * fs  * (gvd*gad);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pINFO) << "f2 : q = " << q2 << ", bar{q} = " << qb2;
    LOG("DISSF", pINFO) << "xf3: q = " << q3 << ", bar{q} = " << qb3;
#endif

    F2val  = q2 + qb2;
    xF3val = 2.0 * (q3-qb3);
  }
  
  // ***  CHARGED CURRENT

  // ***  CHARGED CURRENT
  if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {
      // KA = KV for charm due to its mass
      q    = ( switch_dv * fdv   * ( kV_val_d + kA_val_d )
             + switch_ds * fds   * ( kV_sea_d + kA_sea_d ) ) * fVud2 ;
      q   += ( switch_dv * fdv_c * ( kV_val_d + kV_val_d )
             + switch_ds * fds_c * ( kV_sea_d + kV_sea_d ) ) * fVcd2 ;
      q   +=   switch_s  * fs    * ( kV_sea_s + kA_sea_s )   * fVus2 ;
      q   +=   switch_s  * fs_c  * ( kV_sea_s + kV_sea_s )   * fVcs2;

      qbar  = switch_ubar * fus  * ( kV_sea_u + kA_sea_u ) * fVud2;
      qbar += switch_ubar * fus  * ( kV_sea_u + kA_sea_u ) * fVus2;
      qbar += switch_cbar * fc_c * ( kV_sea_u + kV_sea_u ) * fVcd2;
      qbar += switch_cbar * fc_c * ( kV_sea_u + kV_sea_u ) * fVcs2;
    } else if (is_nubar) {
	    q    = ( switch_uv * fuv  * ( kV_val_u + kA_val_u )
             + switch_us * fus  * ( kV_sea_u + kA_sea_u ) ) * fVud2 ;
	    q   += ( switch_uv * fuv  * ( kV_val_u + kA_val_u )
             + switch_us * fus  * ( kV_sea_u + kA_sea_u ) ) * fVus2 ;
	    q   +=   switch_c  * fc_c * ( kV_sea_u + kV_sea_u )   * fVcd2;
	    q   +=   switch_c  * fc_c * ( kV_sea_u + kV_sea_u )   * fVcs2;

	    qbar  = switch_dbar * fds_c * ( kV_sea_d + kV_sea_d ) * fVcd2;
	    qbar += switch_dbar * fds   * ( kV_sea_d + kA_sea_d ) * fVud2;
	    qbar += switch_sbar * fs    * ( kV_sea_s + kA_sea_s ) * fVus2;
	    qbar += switch_sbar * fs_c  * ( kV_sea_s + kV_sea_s ) * fVcs2;
    } else {
	    return;
    }

    F2val  = (q+qbar);

    if (is_nu) {
      q    = ( switch_dv * fdv   * sqrt( kV_val_d * kA_val_d )
             + switch_ds * fds   * sqrt( kV_sea_d * kA_sea_d ) ) * fVud2;
      q   +=   switch_s  * fs    * sqrt( kV_sea_s * kA_sea_s )   * fVus2;
      q   += ( switch_dv * fdv_c * sqrt( kV_val_d * kV_val_d )
             + switch_ds * fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fVcd2;
      q   +=   switch_s  * fs_c  * sqrt( kV_sea_s * kV_sea_s )   * fVcs2;

      qbar  = switch_ubar * fus  * sqrt( kV_sea_u * kA_sea_u ) * fVud2;
      qbar += switch_ubar * fus  * sqrt( kV_sea_u * kA_sea_u ) * fVus2;
      qbar += switch_cbar * fc_c * sqrt( kV_sea_u * kV_sea_u ) * fVcd2;
      qbar += switch_cbar * fc_c * sqrt( kV_sea_u * kV_sea_u ) * fVcs2;
    }
    else if (is_nubar) {
	    q    = ( switch_uv * fuv * sqrt( kV_val_u * kA_val_u ) + switch_us * fus * sqrt( kV_sea_u * kA_sea_u ) ) * fVud2;
	    q   += ( switch_uv * fuv * sqrt( kV_val_u * kA_val_u ) + switch_us * fus * sqrt( kV_sea_u * kA_sea_u ) ) * fVus2;
	    q   += ( switch_c  * fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fVcd2;
	    q   += ( switch_c  * fc_c * sqrt( kV_sea_u * kV_sea_u ) ) * fVcs2;

	    qbar  = ( switch_dbar * fds_c * sqrt( kV_sea_d * kV_sea_d ) ) * fVcd2;
	    qbar += ( switch_dbar * fds   * sqrt( kV_sea_d * kA_sea_d ) ) * fVud2;
	    qbar += ( switch_sbar * fs    * sqrt( kV_sea_s * kA_sea_s ) ) * fVus2;
	    qbar += ( switch_sbar * fs_c  * sqrt( kV_sea_s * kV_sea_s ) ) * fVcs2;
    } else {
	    return;
    }

    xF3val = 2*(q-qbar);
  }

  // ***  ELECTROMAGNETIC

  if(is_EM) {

    if(!pdg::IsChargedLepton(probe_pdgc)) return;

    double sq23 = TMath::Power(2./3., 2.);
    double sq13 = TMath::Power(1./3., 2.);

    double qu   = sq23 * ( switch_uv   * fuv * kV_val_u + switch_us * fus * kV_sea_u ) ;
    double qd   = sq13 * ( switch_dv   * fdv * kV_val_d + switch_ds * fds * kV_sea_d ) ;
    double qs   = sq13 * ( switch_s    * fs  * kV_sea_s ) ;
    double qbu  = sq23 * ( switch_ubar * fus * kV_sea_u );
    double qbd  = sq13 * ( switch_dbar * fds * kV_sea_d );
    double qbs  = sq13 * ( switch_sbar * fs  * kV_sea_s );

    double q    = qu  + qd  + qs;
    double qbar = qbu + qbd + qbs;

    F2val  = q + qbar;
    xF3val = 0.;

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
//____________________________________________________________________________
double QPMDISStrucFuncBase::Q2(const Interaction * interaction) const
{
// Return Q2 from the kinematics or, if not set, compute it from x,y
// The x might be corrected

  const Kinematics & kinematics = interaction->Kine();

  // if Q2 (or q2) is set then prefer this value
  if (kinematics.KVSet(kKVQ2) || kinematics.KVSet(kKVq2)) {
    double Q2val = kinematics.Q2();
    return Q2val;
  }
  // if Q2 was not set, then compute it from x,y,Ev,Mnucleon
  if (kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->InitState();
    double Mn = init_state.Tgt().HitNucP4Ptr()->M(); // could be off-shell
    //double x     = this->ScalingVar(interaction);       // could be redefined
    double x     = kinematics.x();
    double y     = kinematics.y();
    double Ev    = init_state.ProbeE(kRfHitNucRest);
    double Q2val = 2*Mn*Ev*x*y;
    return Q2val;
  }
  LOG("DISSF", pERROR) << "Could not compute Q2!";
  return 0;
}
//____________________________________________________________________________
double QPMDISStrucFuncBase::ScalingVar(const Interaction* interaction, double Mf) const
{
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  return interaction->Kine().x();
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::KVectorFactors(const Interaction *,
                 double & kuv, double & kdv, double & kus, double & kds, double & ks) const
{
// This is an abstract class: no model-specific correction
// The PDF scaling variables are set to 1
// Override this method to compute model-dependent corrections

  kuv = 1.;
  kdv = 1.;
  kus = 1.;
  kds = 1.;
  ks = 1.;
}

void QPMDISStrucFuncBase::KAxialFactors(const Interaction *,
                 double & kuv, double & kdv, double & kus, double & kds, double & ks) const
{
// This is an abstract class: no model-specific correction
// The PDF scaling variables are set to 1
// Override this method to compute model-dependent corrections

  kuv = 1.;
  kdv = 1.;
  kus = 1.;
  kds = 1.;
  ks = 1.;
}

double QPMDISStrucFuncBase::KCharm(const Interaction * interaction, double Mf) const {
  // Not used by default
  return 1;
}

double QPMDISStrucFuncBase::H(const Interaction * interaction) const {
  // Not used by default
  return 1;
}

//____________________________________________________________________________
double QPMDISStrucFuncBase::R(const Interaction * interaction) const
{
// Computes R ( ~ longitudinal structure function FL = R * 2xF1)
// The scaling variable can be overwritten to include corrections

//   The x used for computing the DIS Nuclear correction factor should be the
//   experimental x, not the rescaled x or off-shell-rest-frame version of x
//   (i.e. selected x).  Since we do not have access to experimental x at this
//   point in the calculation, just use selected x.
  if(fIncludeR) {
    const Kinematics & kine  = interaction->Kine();
    double x  = kine.x();
//    double x  = this->ScalingVar(interaction);
    double Q2val = this->Q2(interaction);
    double Rval  = utils::phys::RWhitlow(x, Q2val);
    return Rval;
  }
  return 0;
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::CalcPDFs(const Interaction * interaction) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  interaction->Print(std::cout);
#endif
  // Clean-up previous calculation
  fPDF  -> Reset();
  fPDFc -> Reset();

  // Get the kinematical variables x,Q2 (could include corrections)
  double x     = this->ScalingVar(interaction);
  double Q2val = this->Q2(interaction);

  // Get the hit nucleon mass (could be off-shell)
  const Target & tgt = interaction->InitState().Tgt();
  double M = tgt.HitNucP4().M();

  // Get the Q2 for which PDFs will be evaluated
  double Q2pdf = TMath::Max(Q2val, fQ2min);

  // Compute PDFs at (x,Q2)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Calculating PDFs @ x = " << x << ", Q2 = " << Q2pdf;
#endif
  fPDF->Calculate(x, Q2pdf);

  // Check whether it is above charm threshold
  bool above_charm = utils::kinematics::IsAboveCharmThreshold(x, Q2val, M, fMc);
  if(above_charm) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pDEBUG)
      << "The event is above the charm threshold (mcharm = " << fMc << ")";
#endif
    if(fCharmOff) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__ 
       LOG("DISSF", pINFO) << "Charm production is turned off";
#endif
    } else {
       // compute the slow rescaling var
       
       double xc = ScalingVar(interaction, fMc);
       if(xc<0 || xc>1) {
          LOG("DISSF", pINFO) << "Unphys. slow rescaling var: xc = " << xc;
       } else {
          // compute PDFs at (xc,Q2)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("DISSF", pDEBUG)
              << "Calculating PDFs @ xc (slow rescaling) = " << x << ", Q2 = " << Q2val;
#endif
          fPDFc->Calculate(xc, Q2pdf);
       }
    }// charm off?
  }//above charm thr?

  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  else {
    LOG("DISSF", pDEBUG)
     << "The event is below the charm threshold (mcharm = " << fMc << ")";
  }
#endif

  // Rules of thumb
  // ---------------------------------------
  // - For W+ exchange use: -1/3|e| quarks and -2/3|e| antiquarks
  // - For W- exchange use:  2/3|e| quarks and  1/3|e| antiquarks
  // - For each qi -> qj transition multiply with the (ij CKM element)^2
  // - Use isospin symmetry to get neutron's u,d from proton's u,d
  //    -- neutron d = proton u
  //    -- neutron u = proton d
  // - Use u = usea + uvalence. Same for d
  // - For s,c use q=qbar
  // - For t,b use q=qbar=0
  double nucScale     = fIncludeNuclMod ? fDISNuclCorr->DISACorrection(interaction) : 1 ;
  fuv   = nucScale * fPDF  -> UpValence();
  fus   = nucScale * fPDF  -> UpSea();
  fdv   = nucScale * fPDF  -> DownValence();
  fds   = nucScale * fPDF  -> DownSea();
  fs    = nucScale * fPDF  -> Strange();
  fc    = nucScale * 0.;
  fuv_c = nucScale * fPDFc -> UpValence();   // will be 0 if < charm threshold
  fus_c = nucScale * fPDFc -> UpSea();       // ...
  fdv_c = nucScale * fPDFc -> DownValence(); // ...
  fds_c = nucScale * fPDFc -> DownSea();     // ...
  fs_c  = nucScale * fPDFc -> Strange();     // ...
  fc_c  = nucScale * fPDFc -> Charm();       // ...
  
  // The above are the proton parton density function. Get the PDFs for the
  // hit nucleon (p or n) by swapping u<->d if necessary

  int nuc_pdgc = tgt.HitNucPdg();
  bool isP = pdg::IsProton  (nuc_pdgc);
  bool isN = pdg::IsNeutron (nuc_pdgc);
  assert(isP  || isN);

  double tmp = 0;
  if (isN) {  // swap u <-> d
    tmp = fuv;   fuv   = fdv;   fdv   = tmp;
    tmp = fus;   fus   = fds;   fds   = tmp;
    tmp = fuv_c; fuv_c = fdv_c; fdv_c = tmp;
    tmp = fus_c; fus_c = fds_c; fds_c = tmp;
  }

}
//____________________________________________________________________________
