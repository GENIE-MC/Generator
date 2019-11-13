//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         This GENIE code was adapted from the neugen3 code co-authored by
         Donna Naples (Pittsburgh U.), Hugh Gallagher (Tufts U), and 
         Costas Andreopoulos (RAL)

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
#include "Physics/BoostedDarkMatter/XSection/QPMDMDISStrucFuncBase.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PhysUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QPMDMDISStrucFuncBase::QPMDMDISStrucFuncBase() :
DISStructureFuncModelI()
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDMDISStrucFuncBase::QPMDMDISStrucFuncBase(string name) :
DISStructureFuncModelI(name)
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDMDISStrucFuncBase::QPMDMDISStrucFuncBase(string name, string config):
DISStructureFuncModelI(name, config)
{
  this->InitPDF();
}
//____________________________________________________________________________
QPMDMDISStrucFuncBase::~QPMDMDISStrucFuncBase()
{
  delete fPDF;
  delete fPDFc;
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::LoadConfig(void)
{
  LOG("DISSF", pDEBUG) << "Loading configuration...";

  //-- pdf
  const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (this->SubAlg("PDF-Set"));
  fPDF  -> SetModel(pdf_model);
  fPDFc -> SetModel(pdf_model);

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

  //-- dark matter couplings
  GetParam( "UpLeftCharge", fQuL );
  GetParam( "UpRightCharge", fQuR );
  GetParam( "CharmLeftCharge", fQcL );
  GetParam( "CharmRightCharge", fQcR );
  GetParam( "DownLeftCharge", fQdL );
  GetParam( "DownRightCharge", fQdR );
  GetParam( "StrangeLeftCharge", fQsL );
  GetParam( "StrangeRightCharge", fQsR );

  LOG("DISSF", pDEBUG) << "Done loading configuration";
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::InitPDF(void)
{
                     // evaluated at:
  fPDF  = new PDF(); //   x = computed (+/-corrections) scaling var, Q2
  fPDFc = new PDF(); //   x = computed charm slow re-scaling var,    Q2
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::Calculate(const Interaction * interaction) const
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
  bool is_dm       = pdg::IsDarkMatter   ( probe_pdgc  );
  bool is_dmb      = pdg::IsAntiDarkMatter   ( probe_pdgc  );
  bool is_dmi      = proc_info.IsDarkMatter();

  if ( !is_dm && !is_dmb    ) return;
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

  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections

  this -> CalcPDFs (interaction);

  //
  // Compute structure functions for the EM, NC and CC cases
  //

  double F2val=0, xF3val=0;

  // ***  DARK MATTER
  if(is_dmi) {

    if(!is_dm && !is_dmb) return;

    double gvu  = 0.5 * (fQuL  + fQuR);
    double gau  = 0.5 * (fQuL  - fQuR);
    double gvc  = 0.5 * (fQcL  + fQcR);
    double gac  = 0.5 * (fQcL  - fQcR);
    double gvd  = 0.5 * (fQdL  + fQdR);
    double gad  = 0.5 * (fQdL  - fQdR);
    double gvs  = 0.5 * (fQsL  + fQsR);
    double gas  = 0.5 * (fQsL  - fQsR);
    double gvu2 = TMath::Power(gvu, 2.);
    double gau2 = TMath::Power(gau, 2.);
    double gvc2 = TMath::Power(gvc, 2.);
    double gac2 = TMath::Power(gac, 2.);
    double gvd2 = TMath::Power(gvd, 2.);
    double gad2 = TMath::Power(gad, 2.);
    double gvs2 = TMath::Power(gvs, 2.);
    double gas2 = TMath::Power(gas, 2.);

    double q2   = 4.0 * ((switch_uv   * fuv + switch_us   * fus) * (gvu2+gau2) + switch_c    * fc  * (gvc2+gac2) + 
			 (switch_dv   * fdv + switch_ds   * fds) * (gvd2+gad2) + switch_s    * fs  * (gvs2+gas2));
    double q3   = 4.0 * ((switch_uv   * fuv + switch_us   * fus) * (2*gvu*gau) + switch_c    * fc  * (2*gvc*gac) + 
			 (switch_dv   * fdv + switch_ds   * fds) * (2*gvd*gad) + switch_s    * fs  * (2*gvs*gas));

    double qb2  = 4.0 * (switch_ubar * fus * (gvu2+gau2) + switch_cbar * fc  * (gvc2+gac2) + 
			 switch_dbar * fds * (gvd2+gad2) + switch_sbar * fs  * (gvs2+gas2));    
    double qb3  = 4.0 * (switch_ubar * fus * (2*gvu*gau) + switch_cbar * fc  * (2*gvc*gac) + 
			 switch_dbar * fds * (2*gvd*gad) + switch_sbar * fs  * (2*gvs*gas));    
 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pINFO) << "f2 : q = " << q2 << ", bar{q} = " << qb2;
    LOG("DISSF", pINFO) << "xf3: q = " << q3 << ", bar{q} = " << qb3;
#endif

    F2val  = q2+qb2;
    xF3val = q3-qb3;
  } 

  double Q2val = this->Q2        (interaction);
  double x     = this->ScalingVar(interaction);
  double f     = this->NuclMod   (interaction); // nuclear modification
  double r     = this->R         (interaction); // R ~ FL

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Nucl. mod   = " << f;
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

    fF3 = f * xF3val/bjx;
    fF2 = f * F2val;
    fF1 = fF2 * 0.5*c/bjx;
    fF5 = fF2/bjx;           // Albright-Jarlskog relation
    fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)
  } 
  else {
    double a = TMath::Power(x,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
    //double a = TMath::Power(x,2.) / Q2val;
    //double c = (1. + 4. * kNucleonMass * a) / (1.+r);

    fF3 = f * xF3val / x;
    fF2 = f * F2val;
    fF1 = fF2 * 0.5 * c / x;
    fF5 = fF2 / x;         // Albright-Jarlskog relation
    fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) 
     << "F1-F5 = " 
     << fF1 << ", " << fF2 << ", " << fF3 << ", " << fF4 << ", " << fF5;
#endif
}
//____________________________________________________________________________
double QPMDMDISStrucFuncBase::Q2(const Interaction * interaction) const
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
double QPMDMDISStrucFuncBase::ScalingVar(const Interaction* interaction) const
{
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  return interaction->Kine().x();
}
//____________________________________________________________________________
void QPMDMDISStrucFuncBase::KFactors(const Interaction *, 
	         double & kuv, double & kdv, double & kus, double & kds) const
{
// This is an abstract class: no model-specific correction
// The PDF scaling variables are set to 1
// Override this method to compute model-dependent corrections

  kuv = 1.;
  kdv = 1.;
  kus = 1.;
  kds = 1.;
}
//____________________________________________________________________________
double QPMDMDISStrucFuncBase::NuclMod(const Interaction * interaction) const
{
// Nuclear modification to Fi
// The scaling variable can be overwritten to include corrections

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return 1.0;
  if( interaction->TestBit(kINoNuclearCorrection) ) return 1.0;

  double f = 1.;
  if(fIncludeNuclMod) {
     const Target & tgt  = interaction->InitState().Tgt();

//   The x used for computing the DIS Nuclear correction factor should be the 
//   experimental x, not the rescaled x or off-shell-rest-frame version of x 
//   (i.e. selected x).  Since we do not have access to experimental x at this 
//   point in the calculation, just use selected x. 
     const Kinematics & kine  = interaction->Kine();
     double x  = kine.x();
     int    A = tgt.A(); 
     f = utils::nuclear::DISNuclFactor(x,A);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DISSF", pDEBUG) << "Nuclear factor for x of " << x << "  = " << f; 
#endif
  }

  return f;
}
//____________________________________________________________________________
double QPMDMDISStrucFuncBase::R(const Interaction * interaction) const
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
void QPMDMDISStrucFuncBase::CalcPDFs(const Interaction * interaction) const
{
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
  bool above_charm = 
           utils::kinematics::IsAboveCharmThreshold(x, Q2val, M, fMc);
  if(above_charm) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pDEBUG) 
      << "The event is above the charm threshold (mcharm = " << fMc << ")";
#endif
    if(fCharmOff) {
       LOG("DISSF", pINFO) << "Charm production is turned off";
    } else {
       // compute the slow rescaling var
       double xc = utils::kinematics::SlowRescalingVar(x, Q2val, M, fMc);    
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
  else {
    LOG("DISSF", pDEBUG) 
     << "The event is below the charm threshold (mcharm = " << fMc << ")";
  }

  // Compute the K factors
  double kval_u = 1.;
  double kval_d = 1.;
  double ksea_u = 1.;
  double ksea_d = 1.;

  this->KFactors(interaction, kval_u, kval_d, ksea_u, ksea_d);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "K-Factors:";
  LOG("DISSF", pDEBUG) << "U: Kval = " << kval_u << ", Ksea = " << ksea_u;
  LOG("DISSF", pDEBUG) << "D: Kval = " << kval_d << ", Ksea = " << ksea_d;
#endif

  // Apply the K factors
  //
  // Always scale d pdfs with d kfactors and u pdfs with u kfactors.
  // Don't swap the applied kfactors for neutrons.
  // Debdatta & Donna noted (Sep.2006) that a similar swap in the neugen
  // implementation was the cause of the difference in nu and nubar F2
  //
  fPDF->ScaleUpValence   (kval_u);
  fPDF->ScaleDownValence (kval_d);
  fPDF->ScaleUpSea       (ksea_u);
  fPDF->ScaleDownSea     (ksea_d);
  fPDF->ScaleStrange     (ksea_d);
  fPDF->ScaleCharm       (ksea_u);
  if(above_charm) {
     fPDFc->ScaleUpValence   (kval_u);
     fPDFc->ScaleDownValence (kval_d);
     fPDFc->ScaleUpSea       (ksea_u);
     fPDFc->ScaleDownSea     (ksea_d);
     fPDFc->ScaleStrange     (ksea_d);
     fPDFc->ScaleCharm       (ksea_u);
  }

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

  fuv   = fPDF  -> UpValence();
  fus   = fPDF  -> UpSea();
  fdv   = fPDF  -> DownValence();
  fds   = fPDF  -> DownSea();
  fs    = fPDF  -> Strange();
  fc    = 0.;
  fuv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  fus_c = fPDFc -> UpSea();       // ...
  fdv_c = fPDFc -> DownValence(); // ...
  fds_c = fPDFc -> DownSea();     // ...
  fs_c  = fPDFc -> Strange();     // ...
  fc_c  = fPDFc -> Charm();       // ...

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
