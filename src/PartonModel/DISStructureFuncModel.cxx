//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 03, 2004

 Adapted from neugen 3.
 Primary authors: D.Naples (Pittsburgh U.), H.Gallagher (Tufts U)

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

  //-- weinberg angle
  double thw = fConfig->GetDoubleDef(
                          "weinberg-angle", gc->GetDouble("WeinbergAngle"));
  fSin2thw = TMath::Power(TMath::Sin(thw), 2);


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

  //-- get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();

  int  nuc_pdgc = tgt.HitNucPdg();
  int  nu_pdgc  = init_state.ProbePdg();
  bool isCC     = proc_info.IsWeakCC();
  bool isNC     = proc_info.IsWeakNC();
  bool isP      = pdg::IsProton       ( nuc_pdgc );
  bool isN      = pdg::IsNeutron      ( nuc_pdgc );
  bool isNu     = pdg::IsNeutrino     ( nu_pdgc  );
  bool isNuBar  = pdg::IsAntiNeutrino ( nu_pdgc  );

  assert(isP  || isN);
  assert(isNu || isNuBar);

  if(tgt.N()==0 && pdg::IsNeutron(nuc_pdgc)) return;
  if(tgt.Z()==0 && pdg::IsProton(nuc_pdgc) ) return;

  if(tgt.HitQrkIsSet()) {
    int  qpdg = tgt.HitQrkPdg();
    bool isu  = pdg::IsUQuark     (qpdg);
    bool isub = pdg::IsAntiUQuark (qpdg);
    bool isd  = pdg::IsDQuark     (qpdg);
    bool isdb = pdg::IsAntiDQuark (qpdg);
    bool iss  = pdg::IsSQuark     (qpdg);
    bool issb = pdg::IsAntiSQuark (qpdg);

    bool isNuCC    = isNu    && isCC;
    bool isNuBarCC = isNuBar && isCC;
    if(isNuCC    && isu ) return;
    if(isNuCC    && isdb) return;
    if(isNuCC    && issb) return;
    if(isNuBarCC && isub) return;
    if(isNuBarCC && isd ) return;
    if(isNuBarCC && iss ) return;
  }
   
  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the rescaling variable and the
  // K factors are applied
  //
  this -> CalcPDFs (interaction);

  // Compute structure functions
  //
  double F2=0, xF3=0;

  // ***  NEUTRAL CURRENT
  //
  if(isNC) {
    double GL   = (isNu) ? ( 0.5 - (2./3.)*fSin2thw) : (     - (2./3.)*fSin2thw); // clu
    double GR   = (isNu) ? (     - (2./3.)*fSin2thw) : ( 0.5 - (2./3.)*fSin2thw); // cru
    double GLp  = (isNu) ? (-0.5 + (1./3.)*fSin2thw) : (       (1./3.)*fSin2thw); // cld
    double GRp  = (isNu) ? (       (1./3.)*fSin2thw) : (-0.5 + (1./3.)*fSin2thw); // crd
    double gvu  = GL  + GR;
    double gau  = GL  - GR;
    double gvd  = GLp + GRp;
    double gad  = GLp - GRp;
    double gvu2 = TMath::Power(gvu, 2.);
    double gau2 = TMath::Power(gau, 2.);
    double gvd2 = TMath::Power(gvd, 2.);
    double gad2 = TMath::Power(gad, 2.);
    double fc   = 0;
    double q2   = (fu+fc)  * (gvu2+gau2) + (fd+fs)  * (gvd2+gad2);
    double q3   = (fu+fc)  * (2*gvu*gau) + (fd+fs)  * (2*gvd*gad);
    double qb2  = (fus+fc) * (gvu2+gau2) + (fds+fs) * (gvd2+gad2);    
    double qb3  = (fus+fc) * (2*gvu*gau) + (fds+fs) * (2*gvd*gad);    
 
    if(tgt.HitQrkIsSet()) {
       // \bar{q} contributions are computed from q(sea) pdfs.
       // Explicity zero q contributions if we have a hit anti-quark.
       // Vice-versa for hit-quarks.
       // One has to do that despite the pdf zeroing in CalcPDFs() as the non-zero
       // q(sea) pdf would give non-zero contributions to S/F from both q and \bar{q}.
       int qpdg = tgt.HitQrkPdg();
       q2  *= ( pdg::IsQuark(qpdg)     ? 1. : 0. );
       q3  *= ( pdg::IsQuark(qpdg)     ? 1. : 0. );
       qb2 *= ( pdg::IsAntiQuark(qpdg) ? 1. : 0. );
       qb3 *= ( pdg::IsAntiQuark(qpdg) ? 1. : 0. );
    }

    LOG("DISSF", pINFO) << "f2 : q = " << q2 << ", bar{q} = " << qb2;
    LOG("DISSF", pINFO) << "xf3: q = " << q3 << ", bar{q} = " << qb3;

    F2  = q2+qb2;
    xF3 = q3-qb3;
  } 

  // ***  CHARGED CURRENT
  //
  if(isCC) {
    double q=0, qbar=0;
    if (isNu) {
      q    = (fd  * fVud2) + (fs  * fVus2) + (fd_c * fVcd2) + (fs_c * fVcs2);
      qbar = (fus * fVud2) + (fus * fVus2) + (fc_c * fVcd2) + (fc_c * fVcs2);
    }
    else if (isNuBar) {
      q    = (fu    * fVud2) + (fu  * fVus2) + (fc_c * fVcd2) + (fc_c * fVcs2);
      qbar = (fds_c * fVcd2) + (fds * fVud2) + (fs   * fVus2) + (fs_c * fVcs2);
    }
    LOG("DISSF", pINFO) << "Q(x,Q2) = " << q << ", Qbar(x,Q2) = " << qbar;
    F2  = 2*(q+qbar);
    xF3 = 2*(q-qbar);
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

  fF3 = f * xF3/x;
  fF2 = f * F2;
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

  if( interaction->TestBit(kIAssumeFreeNucleon)   ) return 1.0;
  if( interaction->TestBit(kINoNuclearCorrection) ) return 1.0;

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

  // Take the sums of valence & sea pdfs
  fu    = fuv   + fus;
  fd    = fdv   + fds;
  fu_c  = fuv_c + fus_c;
  fd_c  = fdv_c + fds_c;

  // The above can be used to compute vN->lX cross sections taking into account 
  // the contribution from all quarks. Check whether a struck quark has been set: 
  // in this case the conributions from all other quarks will be set to 0 as as 
  // this algorithm can be used from computing vq->lq cross sections as well.

  if(tgt.HitQrkIsSet()) {
    int  qpdg = tgt.HitQrkPdg();
    bool sea  = tgt.HitSeaQrk();

    bool isu  = pdg::IsUQuark     (qpdg);
    bool isub = pdg::IsAntiUQuark (qpdg);
    bool isd  = pdg::IsDQuark     (qpdg);
    bool isdb = pdg::IsAntiDQuark (qpdg);
    bool iss  = pdg::IsSQuark     (qpdg);
    bool issb = pdg::IsAntiSQuark (qpdg);

    fuv   = ( isu        && !sea) ? fuv   : 0.;
    fus   = ((isu||isub) &&  sea) ? fus   : 0.; 
    fdv   = ( isd        && !sea) ? fdv   : 0.;
    fds   = ((isd||isdb) &&  sea) ? fds   : 0.;
    fs    = ((iss||issb) &&  sea) ? fs    : 0.;
    fuv_c = ( isu        && !sea) ? fuv_c : 0.;
    fus_c = ((isu||isub) &&  sea) ? fus_c : 0.;
    fdv_c = ( isd        && !sea) ? fdv_c : 0.;
    fds_c = ((isd||isdb) &&  sea) ? fds_c : 0.;
    fs_c  = ((iss||issb) &&  sea) ? fs_c  : 0.;
    fc_c  = 0;

    fu    = fuv   + fus;
    fd    = fdv   + fds;
    fu_c  = fuv_c + fus_c;
    fd_c  = fdv_c + fds_c;
  }
/*
  LOG("DISSF", pDEBUG) << "u(v) = " << fuv;
  LOG("DISSF", pDEBUG) << "u(s) = " << fus;
  LOG("DISSF", pDEBUG) << "d(v) = " << fdv;
  LOG("DISSF", pDEBUG) << "d(s) = " << fds;
  LOG("DISSF", pDEBUG) << "s    = " << fs;
*/
}
//____________________________________________________________________________
