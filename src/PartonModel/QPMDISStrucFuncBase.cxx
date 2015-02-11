//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 Adapted from neugen 3.
 Primary authors: D.Naples (Pittsburgh U.), H.Gallagher (Tufts U), CA

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 08, 2009 - CA
   Added structure function calculations for charged leptons.
   Restructured to make switching off quark contributions and going from 
   lN->l'X mode to lq->l'q' mode more transparant
 @ Oct 09, 2009 - CA
   Renamed to QPMDISStrucFuncBase from DISStructureFuncModel.
 @ Aug 12, 2014 - HG 
   Fix a problem identified by Brian Tice (Minerva)
   The nuclear modification to the pdf should be calculated in terms 
   of the experimental x, not the rescaled x.  The same goes for R(x,Q2).
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Messenger/Messenger.h"
#include "PartonModel/QPMDISStrucFuncBase.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"
#include "Utils/PhysUtils.h"

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

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- pdf
  const PDFModelI * pdf_model =
         dynamic_cast<const PDFModelI *> (this->SubAlg("PDF-Set"));
  fPDF  -> SetModel(pdf_model);
  fPDFc -> SetModel(pdf_model);

  //-- get CKM elements
  fVcd  = fConfig->GetDoubleDef("CKM-Vcd", gc->GetDouble("CKM-Vcd"));
  fVcs  = fConfig->GetDoubleDef("CKM-Vcs", gc->GetDouble("CKM-Vcs"));
  fVud  = fConfig->GetDoubleDef("CKM-Vud", gc->GetDouble("CKM-Vud"));
  fVus  = fConfig->GetDoubleDef("CKM-Vus", gc->GetDouble("CKM-Vus"));

  fVcd2 = TMath::Power( fVcd, 2 );
  fVcs2 = TMath::Power( fVcs, 2 );
  fVud2 = TMath::Power( fVud, 2 );
  fVus2 = TMath::Power( fVus, 2 );

  //-- charm mass
  fMc = fConfig->GetDoubleDef("Charm-Mass", gc->GetDouble("Charm-Mass"));

  //-- min Q2 for PDF evaluation
  fQ2min = fConfig->GetDoubleDef("PDF-Q2min", gc->GetDouble("PDF-Q2min"));

  //-- include R (~FL)?
  fIncludeR = fConfig->GetBoolDef(
                           "IncludeR", gc->GetBool("DISSF-IncludeR"));

  //-- include nuclear factor (shadowing / anti-shadowing / ...)
  fIncludeNuclMod = fConfig->GetBoolDef(
               "IncludeNuclMod", gc->GetBool("DISSF-IncludeNuclMod"));

  //-- turn charm production off?
  fCharmOff  = fConfig->GetBoolDef("Charm-Prod-Off", false);

  //-- weinberg angle
  double thw = fConfig->GetDoubleDef(
                          "weinberg-angle", gc->GetDouble("WeinbergAngle"));
  fSin2thw = TMath::Power(TMath::Sin(thw), 2);

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
  bool is_CC       = proc_info.IsWeakCC();
  bool is_NC       = proc_info.IsWeakNC();
  bool is_EM       = proc_info.IsEM();

  if ( !is_lepton           ) return;
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

  this -> CalcPDFs (interaction);

  //
  // Compute structure functions for the EM, NC and CC cases
  //

  double F2=0, xF3=0;

  // ***  NEUTRAL CURRENT

  if(is_NC) {

    if(!is_nu && !is_nubar) return;

    double GL   = (is_nu) ? ( 0.5 - (2./3.)*fSin2thw) : (     - (2./3.)*fSin2thw); // clu
    double GR   = (is_nu) ? (     - (2./3.)*fSin2thw) : ( 0.5 - (2./3.)*fSin2thw); // cru
    double GLp  = (is_nu) ? (-0.5 + (1./3.)*fSin2thw) : (       (1./3.)*fSin2thw); // cld
    double GRp  = (is_nu) ? (       (1./3.)*fSin2thw) : (-0.5 + (1./3.)*fSin2thw); // crd
    double gvu  = GL  + GR;
    double gau  = GL  - GR;
    double gvd  = GLp + GRp;
    double gad  = GLp - GRp;
    double gvu2 = TMath::Power(gvu, 2.);
    double gau2 = TMath::Power(gau, 2.);
    double gvd2 = TMath::Power(gvd, 2.);
    double gad2 = TMath::Power(gad, 2.);

    double q2   = (switch_uv   * fuv + switch_us   * fus + switch_c    * fc)  * (gvu2+gau2) + 
                  (switch_dv   * fdv + switch_ds   * fds + switch_s    * fs)  * (gvd2+gad2);
    double q3   = (switch_uv   * fuv + switch_us   * fus + switch_c    * fc)  * (2*gvu*gau) + 
                  (switch_dv   * fdv + switch_ds   * fds + switch_s    * fs)  * (2*gvd*gad);

    double qb2  = (switch_ubar * fus + switch_cbar * fc)  * (gvu2+gau2) + 
                  (switch_dbar * fds + switch_sbar * fs)  * (gvd2+gad2);    
    double qb3  = (switch_ubar * fus + switch_cbar * fc)  * (2*gvu*gau) + 
                  (switch_dbar * fds + switch_sbar * fs)  * (2*gvd*gad);    
 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pINFO) << "f2 : q = " << q2 << ", bar{q} = " << qb2;
    LOG("DISSF", pINFO) << "xf3: q = " << q3 << ", bar{q} = " << qb3;
#endif

    F2  = q2+qb2;
    xF3 = q3-qb3;
  } 

  // ***  CHARGED CURRENT

  if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {
      q    = ( switch_dv * fdv   + switch_ds * fds   ) * fVud2 + 
             ( switch_s  * fs                        ) * fVus2 + 
             ( switch_dv * fdv_c + switch_ds * fds_c ) * fVcd2 + 
             ( switch_s  * fs_c                      ) * fVcs2;

      qbar = ( switch_ubar * fus  ) * fVud2 + 
             ( switch_ubar * fus  ) * fVus2 + 
             ( switch_cbar * fc_c ) * fVcd2 + 
             ( switch_cbar * fc_c ) * fVcs2;
    }
    else 
    if (is_nubar) {
      q    = ( switch_uv * fuv + switch_us * fus    ) * fVud2 + 
             ( switch_uv * fuv + switch_us * fus    ) * fVus2 + 
             ( switch_c  * fc_c                     ) * fVcd2 + 
             ( switch_c  * fc_c                     ) * fVcs2;

      qbar = ( switch_dbar * fds_c ) * fVcd2 + 
             ( switch_dbar * fds   ) * fVud2 + 
             ( switch_sbar * fs    ) * fVus2 + 
             ( switch_sbar * fs_c  ) * fVcs2;
    }
    else {
      return;
    }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pINFO) << "Q(x,Q2) = " << q << ", Qbar(x,Q2) = " << qbar;
#endif
    F2  = 2*(q+qbar);
    xF3 = 2*(q-qbar);
  }

  // ***  ELECTROMAGNETIC

  if(is_EM) {

    if(!pdg::IsChargedLepton(probe_pdgc)) return;

    double sq23 = TMath::Power(2./3., 2.);
    double sq13 = TMath::Power(1./3., 2.);

    double qu   = sq23 * ( switch_uv   * fuv + switch_us * fus );
    double qd   = sq13 * ( switch_dv   * fdv + switch_ds * fds );
    double qs   = sq13 * ( switch_s    * fs  );
    double qbu  = sq23 * ( switch_ubar * fus );
    double qbd  = sq13 * ( switch_dbar * fds );
    double qbs  = sq13 * ( switch_sbar * fs  );

    double q    = qu  + qd  + qs;
    double qbar = qbu + qbd + qbs;

    F2  = q + qbar;;
    xF3 = 0.;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pINFO) << "Q(x,Q2) = " << q << ", Qbar(x,Q2) = " << qbar;
#endif

  }

  double Q2 = this->Q2        (interaction);
  double x  = this->ScalingVar(interaction);
  double f  = this->NuclMod   (interaction); // nuclear modification
  double r  = this->R         (interaction); // R ~ FL

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Nucl. mod   = " << f;
  LOG("DISSF", pDEBUG) << "R(=FL/2xF1) = " << r;
#endif

  double a = TMath::Power(x,2.) / TMath::Max(Q2, 0.8);
  double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
//double a = TMath::Power(x,2.) / Q2;
//double c = (1. + 4. * kNucleonMass * a) / (1.+r);

  fF3 = f * xF3/x;
  fF2 = f * F2;
  fF1 = fF2 * 0.5*c/x;
  fF5 = fF2/x;           // Albright-Jarlskog relation
  fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)

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
    double Q2 = kinematics.Q2();
    return Q2;
  }
  // if Q2 was not set, then compute it from x,y,Ev,Mnucleon
  if (kinematics.KVSet(kKVy)) {
    const InitialState & init_state = interaction->InitState();
    double Mn = init_state.Tgt().HitNucP4Ptr()->M(); // could be off-shell
    //double x  = this->ScalingVar(interaction);       // could be redefined
    double x  = kinematics.x();
    double y  = kinematics.y();
    double Ev = init_state.ProbeE(kRfHitNucRest);
    double Q2 = 2*Mn*Ev*x*y;
    return Q2;
  }
  LOG("DISSF", pERROR) << "Could not compute Q2!";
  return 0;
}
//____________________________________________________________________________
double QPMDISStrucFuncBase::ScalingVar(const Interaction* interaction) const
{
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  return interaction->Kine().x();
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::KFactors(const Interaction *, 
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
double QPMDISStrucFuncBase::NuclMod(const Interaction * interaction) const
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
    double Q2 = this->Q2(interaction);
    double R = utils::phys::RWhitlow(x, Q2);
    return R;
  }
  return 0;
}
//____________________________________________________________________________
void QPMDISStrucFuncBase::CalcPDFs(const Interaction * interaction) const
{
  // Clean-up previous calculation
  fPDF  -> Reset();
  fPDFc -> Reset();

  // Get the kinematical variables x,Q2 (could include corrections)
  double x  = this->ScalingVar(interaction);
  double Q2 = this->Q2(interaction);

  // Get the hit nucleon mass (could be off-shell)
  const Target & tgt = interaction->InitState().Tgt();
  double M = tgt.HitNucP4().M(); 

  // Get the Q2 for which PDFs will be evaluated
  double Q2pdf = TMath::Max(Q2, fQ2min);

  // Compute PDFs at (x,Q2)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISSF", pDEBUG) << "Calculating PDFs @ x = " << x << ", Q2 = " << Q2pdf;
#endif
  fPDF->Calculate(x, Q2pdf);

  // Check whether it is above charm threshold
  bool above_charm = 
           utils::kinematics::IsAboveCharmThreshold(x,Q2,M,fMc);
  if(above_charm) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("DISSF", pDEBUG) 
      << "The event is above the charm threshold (mcharm = " << fMc << ")";
#endif
    if(fCharmOff) {
       LOG("DISSF", pINFO) << "Charm production is turned off";
    } else {
       // compute the slow rescaling var
       double xc = utils::kinematics::SlowRescalingVar(x,Q2,M,fMc);    
       if(xc<0 || xc>1) {
          LOG("DISSF", pINFO) << "Unphys. slow rescaling var: xc = " << xc;
       } else {
          // compute PDFs at (xc,Q2)
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("DISSF", pDEBUG) 
              << "Calculating PDFs @ xc (slow rescaling) = " << x << ", Q2 = " << Q2;
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
