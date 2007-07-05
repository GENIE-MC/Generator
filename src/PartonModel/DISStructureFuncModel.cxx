//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModel.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel() :
DISStructureFuncModelI()
{
  this->InitPDF();
}
//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel(string name) :
DISStructureFuncModelI(name)
{
  this->InitPDF();
}
//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel(string name, string config):
DISStructureFuncModelI(name, config)
{
  this->InitPDF();
}
//____________________________________________________________________________
DISStructureFuncModel::~DISStructureFuncModel()
{
  delete fPDF;
  delete fPDFc;
}
//____________________________________________________________________________
void DISStructureFuncModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISStructureFuncModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISStructureFuncModel::LoadConfig(void)
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

  LOG("DISSF", pDEBUG) << "Done loading configuration";
}
//____________________________________________________________________________
void DISStructureFuncModel::InitPDF(void)
{
                     // at each calculation are evaluated at:
  fPDF  = new PDF(); //   x = computed (+/-corrections) scaling var, Q2
  fPDFc = new PDF(); //   x = computed charm slow re-scaling var,    Q2
}
//____________________________________________________________________________
void DISStructureFuncModel::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  // Some tests in case it is called directly & not from an XSecAlgorithmI
  // which has validated the input interaction
  if(!interaction->TestBit(kISkipProcessChk)) {
    const Target & tgt = interaction->InitState().Tgt();
    if(!tgt.HitNucIsSet()) return;
    int nuc = tgt.HitNucPdg();
    if(tgt.N()==0 && pdg::IsNeutron(nuc)) return;
    if(tgt.Z()==0 && pdg::IsProton(nuc) ) return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the rescaling variable and the
  // K factors are applied

  double q    = 0;
  double qbar = 0;

  this -> CalcPDFs (interaction);
  this -> QQBar    (interaction,q,qbar);

  if(q<0 || qbar<0) {
     LOG("DISSF", pERROR) << "Negative q and/or q{bar}! Can not compute SFs";
     return;
  }
  double Q2 = this->Q2        (interaction);
  double x  = this->ScalingVar(interaction);
  double f  = this->NuclMod   (interaction); // nuclear modification
  double r  = this->R         (interaction); // R ~ FL

  LOG("DISSF", pDEBUG) << "Nucl. mod   = " << f;
  LOG("DISSF", pDEBUG) << "R(=FL/2xF1) = " << r;

  double a = TMath::Power(x,2.) / TMath::Max(Q2, 0.8);
  double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);

  //double a = TMath::Power(x,2.) / Q2;
  //double c = (1. + 4. * kNucleonMass * a) / (1.+r);

  fF3 = f * 2*(q-qbar)/x;
  fF2 = f * 2*(q+qbar);
  fF1 = fF2 * 0.5*c/x;
  fF5 = fF2/x;           // Albright-Jarlskog relation
  fF4 = 0.;              // Nucl.Phys.B 84, 467 (1975)

  LOG("DISSF", pDEBUG) 
     << "F1-F5 = " 
     << fF1 << ", " << fF2 << ", " << fF3 << ", " << fF4 << ", " << fF5;
}
//____________________________________________________________________________
double DISStructureFuncModel::Q2(const Interaction * interaction) const
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
double DISStructureFuncModel::ScalingVar(const Interaction* interaction) const
{
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  return interaction->Kine().x();
}
//____________________________________________________________________________
void DISStructureFuncModel::KFactors(const Interaction *, 
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
double DISStructureFuncModel::NuclMod(const Interaction * interaction) const
{
// Nuclear modification to Fi
// The scaling variable can be overwritten to include corrections

  // if requested switch off nuclear corrections even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return 1.0;

  double f = 1.;
  if(fIncludeNuclMod) {
     const Target & tgt  = interaction->InitState().Tgt();
     double x = this->ScalingVar(interaction);
     int    A = tgt.A();
     f = utils::nuclear::DISNuclFactor(x,A);
  }
  return f;
}
//____________________________________________________________________________
double DISStructureFuncModel::R(const Interaction * interaction) const
{
// Computes R ( ~ longitudinal structure function FL = R * 2xF1)
// The scaling variable can be overwritten to include corrections

  if(fIncludeR) {
    double x  = this->ScalingVar(interaction);
    double Q2 = this->Q2(interaction);
    double R = utils::nuclear::RModelMod(x, Q2);
    return R;
  }
  return 0;
}
//____________________________________________________________________________
void DISStructureFuncModel::CalcPDFs(const Interaction * interaction) const
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
  LOG("DISSF", pDEBUG) << "Calculating PDFs @ x = " << x << ", Q2 = " << Q2pdf;
  fPDF->Calculate(x, Q2pdf);

  // Check whether it is above charm threshold
  bool above_charm = 
           utils::kinematics::IsAboveCharmThreshold(x,Q2,M,fMc);
  if(above_charm) {
    LOG("DISSF", pDEBUG) 
      << "The event is above the charm threshold (mcharm = " << fMc << ")";

    if(fCharmOff) {
       LOG("DISSF", pINFO) << "Charm production is turned off";
    } else {
       // compute the slow rescaling var
       double xc = utils::kinematics::SlowRescalingVar(x,Q2,M,fMc);    
       if(xc<0 || xc>1) {
          LOG("DISSF", pINFO) << "Unphys. slow rescaling var: xc = " << xc;
       } else {
          // compute PDFs at (xc,Q2)
          LOG("DISSF", pDEBUG) 
              << "Calculating PDFs @ xc (slow rescaling) = " 
              << x << ", Q2 = " << Q2;
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

  LOG("DISSF", pDEBUG) << "K-Factors:";
  LOG("DISSF", pDEBUG) << "U: Kval = " << kval_u << ", Ksea = " << ksea_u;
  LOG("DISSF", pDEBUG) << "D: Kval = " << kval_d << ", Ksea = " << ksea_d;

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
}
//____________________________________________________________________________
void DISStructureFuncModel::QQBar(
             const Interaction * interaction, double & q, double & qbar) const
{
  q     = 0;
  qbar  = 0;

  double tmp = 0;

  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();

  int nuc_pdgc = target.HitNucPdg();
  int nu_pdgc  = init_state.ProbePdg();

  bool isP     = pdg::IsProton       ( nuc_pdgc );
  bool isN     = pdg::IsNeutron      ( nuc_pdgc );
  bool isNu    = pdg::IsNeutrino     ( nu_pdgc  );
  bool isNuBar = pdg::IsAntiNeutrino ( nu_pdgc  );

  assert(isP  || isN);
  assert(isNu || isNuBar);

  // Rules of thumb for computing Q and QBar
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

  // Get PDFs [should have been computed by calling in CalcPDFs() first]

  double uv   = fPDF  -> UpValence();
  double us   = fPDF  -> UpSea();
  double dv   = fPDF  -> DownValence();
  double ds   = fPDF  -> DownSea();
  double s    = fPDF  -> Strange();
  double uv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  double us_c = fPDFc -> UpSea();       // ...
  double dv_c = fPDFc -> DownValence(); // ...
  double ds_c = fPDFc -> DownSea();     // ...
  double s_c  = fPDFc -> Strange();     // ...
  double c_c  = fPDFc -> Charm();       // ...

  // The above are the proton parton density function. Get the PDFs for the 
  // hit nucleon (p or n) by swapping u<->d if necessary

  if (isN) {  // swap u <-> d
    tmp = uv;   uv   = dv;   dv   = tmp;
    tmp = us;   us   = ds;   ds   = tmp;
    tmp = uv_c; uv_c = dv_c; dv_c = tmp;
    tmp = us_c; us_c = ds_c; ds_c = tmp;
  }

  // Take the sums of valence & sea pdfs
  double u    = uv   + us;
  double d    = dv   + ds;
  double u_c  = uv_c + us_c;
  double d_c  = dv_c + ds_c;

  // The above can be used to compute vN->lX cross sections taking into account 
  // the contribution from all quarks. Check whether a struck quark has been set: 
  // in this case the conributions from all other quarks will be set to 0 as as 
  // this algorithm can be used from computing vq->lq cross sections as well.

  if(target.HitQrkIsSet()) {
    int  qpdg = target.HitQrkPdg();
    bool sea  = target.HitSeaQrk();

    bool isu  = pdg::IsUQuark     (qpdg);
    bool isub = pdg::IsAntiUQuark (qpdg);
    bool isd  = pdg::IsDQuark     (qpdg);
    bool isdb = pdg::IsAntiDQuark (qpdg);
    bool iss  = pdg::IsSQuark     (qpdg);
    bool issb = pdg::IsAntiSQuark (qpdg);

    const ProcessInfo & proc_info = interaction->ProcInfo();
    bool isNuCC    = isNu    && proc_info.IsWeakCC();
    bool isNuBarCC = isNuBar && proc_info.IsWeakCC();
    if(isNuCC    && isu ) return;
    if(isNuCC    && isdb) return;
    if(isNuCC    && issb) return;
    if(isNuBarCC && isub) return;
    if(isNuBarCC && isd ) return;
    if(isNuBarCC && iss ) return;

    uv   = ( isu        && !sea) ? uv   : 0.;
    us   = ((isu||isub) &&  sea) ? us   : 0.; 
    dv   = ( isd        && !sea) ? dv   : 0.;
    ds   = ((isd||isdb) &&  sea) ? ds   : 0.;
    s    = ((iss||issb) &&  sea) ? s    : 0.;
    uv_c = ( isu        && !sea) ? uv_c : 0.;
    us_c = ((isu||isub) &&  sea) ? us_c : 0.;
    dv_c = ( isd        && !sea) ? dv_c : 0.;
    ds_c = ((isd||isdb) &&  sea) ? ds_c : 0.;
    s_c  = ((iss||issb) &&  sea) ? s_c  : 0.;
    c_c  = 0;

    u    = uv   + us;
    d    = dv   + ds;
    u_c  = uv_c + us_c;
    d_c  = dv_c + ds_c;
  }

  // Compute q, qbar for the contributing quarks
  if (isNu) {
    q    = (d  * fVud2) + (s  * fVus2) + (d_c * fVcd2) + (s_c * fVcs2);
    qbar = (us * fVud2) + (us * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
  }
  else if (isNuBar) {
    q    = (u    * fVud2) + (u  * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
    qbar = (ds_c * fVcd2) + (ds * fVud2) + (s   * fVus2) + (s_c * fVcs2);
  }
  else {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }

  LOG("DISSF", pDEBUG) << "Q(x,Q2) = " << q << ", Qbar(x,Q2) = " << qbar;
}
//____________________________________________________________________________
