//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
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
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- pdf
  const PDFModelI * pdf_model =
                  dynamic_cast<const PDFModelI *>
                             (this->SubAlg("pdf-alg-name", "pdf-param-set"));
  fPDF  -> SetModel(pdf_model);
  fPDFc -> SetModel(pdf_model);

  //-- get CKM elements
  fVcd  = fConfig->GetDoubleDef("Vcd", gc->GetDouble("CKM-Vcd"));
  fVcs  = fConfig->GetDoubleDef("Vcs", gc->GetDouble("CKM-Vcs"));
  fVud  = fConfig->GetDoubleDef("Vud", gc->GetDouble("CKM-Vud"));
  fVus  = fConfig->GetDoubleDef("Vus", gc->GetDouble("CKM-Vus"));

  fVcd2 = TMath::Power( fVcd, 2 );
  fVcs2 = TMath::Power( fVcs, 2 );
  fVud2 = TMath::Power( fVud, 2 );
  fVus2 = TMath::Power( fVus, 2 );

  //-- charm mass
  fMc = fConfig->GetDoubleDef("c-quark-mass", gc->GetDouble("Charm-Mass"));

  //-- min Q2 for PDF evaluation
  fQ2min = fConfig->GetDoubleDef("Q2min", gc->GetDouble("PDF-Q2min"));


  fIncludeFL      = true;
  fIncludeNuclMod = true;
  fCorrectF3      = true;
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

  double Q2 = this->Q2(interaction);

  const Kinematics & kine = interaction->GetKinematics();
  double x = kine.x();
  if(x<=0.) {
     LOG("DISSF", pERROR)
             << "scaling variable x = " << x << " < 0. Can not compute SFs";
     return;
  }

  double f = this->NuclMod (interaction); // nuclear modification
  double r = this->R       (interaction); // R ~ FL

  LOG("DISSF", pDEBUG) << "Nucl. mod   = " << f;
  LOG("DISSF", pDEBUG) << "R(=FL/2xF1) = " << r;

  fF3 = f * 2*(q-qbar)/x;
  fF2 = f * 2*(q+qbar);

  if(fCorrectF3) {
    double a = TMath::Power(x,2.) / TMath::Max(Q2, 0.8);
    double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);
    fF3 = fF3 * c;
    fF1 = fF2 * 0.5*c/x;
  }
  else {
    double a = TMath::Power(x,2.) / Q2;
    double c = (1. + 4. * kNucleonMass * a) / (1.+r);
    fF1 = fF2 * 0.5*c/x;
  }

  fF5 = fF2/x; // Albright-Jarlskog relations
  fF4 = 0.;    // Nucl.Phys.B 84, 467 (1975)
}
//____________________________________________________________________________
double DISStructureFuncModel::Q2(const Interaction * interaction) const
{
// Return Q2 from the kinematics or, if not set, compute it from x,y

  return utils::kinematics::CalcQ2(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModel::ScalingVar(
                                        const Interaction * interaction) const
{
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  return interaction->GetKinematics().x();
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

  // if requested switch off nuclear corrections even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return 1.0;

  double f = 1.;

  if(fIncludeNuclMod) {
     const Kinematics & kine = interaction->GetKinematics();
     const Target &     tgt  = interaction->GetInitialState().GetTarget();
     double x = kine.x();
     int    A = tgt.A();
     f = utils::nuclear::DISNuclFactor(x,A);
  }
  return f;
}
//____________________________________________________________________________
double DISStructureFuncModel::R(const Interaction * interaction) const
{
// Computes R ( ~ longitudinal structure function FL = R * 2xF1)

  double R = 0;

  if(fIncludeFL) {
    double x  = interaction->GetKinematics().x();
    double Q2 = this->Q2(interaction);//Q2 from kinematics or compute from x,y
    R = utils::nuclear::RModelMod(x, Q2);
  }
  return R;
}
//____________________________________________________________________________
void DISStructureFuncModel::CalcPDFs(const Interaction * interaction) const
{
  //-- Clean-up previous calculation
  fPDF  -> Reset();
  fPDFc -> Reset();

  //-- Get the kinematical variables (x,Q2)
  double x  = this->ScalingVar(interaction);
  double Q2 = this->Q2(interaction);

  //-- apply Q2 cut
  Q2 = TMath::Max(Q2, fQ2min);

  //-- compute PDFs at (x,Q2)
  fPDF->Calculate(x, Q2);

  //-- check whether it is above charm threshold
  bool isAbvCh = utils::kinematics::IsAboveCharmThreshold(interaction, fMc);
  if(isAbvCh) {
    // compute PDFs at (xi,Q2)
    double xc = utils::kinematics::SlowRescalingVar(interaction, fMc);    
    if(xc>0 && xc<1) fPDFc->Calculate(xc, Q2);
  }

  //-- compute the K factors
  double kval_u = 1.;
  double kval_d = 1.;
  double ksea_u = 1.;
  double ksea_d = 1.;

  this->KFactors(interaction, kval_u, kval_d, ksea_u, ksea_d);

  LOG("DISSF", pDEBUG) << "U: Kval = " << kval_u << ", Ksea = " << ksea_u;
  LOG("DISSF", pDEBUG) << "D: Kval = " << kval_d << ", Ksea = " << ksea_d;

  //-- apply the K factors
  fPDF->ScaleUpValence   (kval_u);
  fPDF->ScaleDownValence (kval_d);
  fPDF->ScaleUpSea       (ksea_u);
  fPDF->ScaleDownSea     (ksea_d);
  if(isAbvCh) {
    fPDFc->ScaleUpValence   (kval_u);
    fPDFc->ScaleDownValence (kval_d);
    fPDFc->ScaleUpSea       (ksea_u);
    fPDFc->ScaleDownSea     (ksea_d);
  }
}
//____________________________________________________________________________
void DISStructureFuncModel::QQBar(
             const Interaction * interaction, double & q, double & qbar) const
{
  q     = -1;
  qbar  = -1;

  const InitialState & init_state = interaction->GetInitialState();

  int nuc_pdgc = init_state.GetTarget().StruckNucleonPDGCode();
  int nu_pdgc  = init_state.GetProbePDGCode();
  bool isP     = pdg::IsProton       ( nuc_pdgc );
  bool isN     = pdg::IsNeutron      ( nuc_pdgc );
  bool isNu    = pdg::IsNeutrino     ( nu_pdgc  );
  bool isNuBar = pdg::IsAntiNeutrino ( nu_pdgc  );

  //-- get PDFs [should have been computed by calling in CalcPDFs() first]

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
  double u    = uv   + us;
  double d    = dv   + ds;
  double u_c  = uv_c + us_c;
  double d_c  = dv_c + ds_c;

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

  if (isP && isNu)
  {
    q    = (d  * fVud2) + (s  * fVus2) + (d_c * fVcd2) + (s_c * fVcs2);
    qbar = (us * fVud2) + (us * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
  }
  else if (isP && isNuBar)
  {
    q    = (u    * fVud2) + (u  * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
    qbar = (ds_c * fVcd2) + (ds * fVud2) + (s   * fVus2) + (s_c * fVcs2);
  }
  else if (isN && isNu)
  {
    q    = (u  * fVud2) + (s  * fVus2) + (u_c * fVcd2) + (s_c * fVcs2);
    qbar = (ds * fVud2) + (ds * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
  }
  else if (isN && isNuBar)
  {
     q   = (d    * fVud2) + (d  * fVus2) + (c_c * fVcd2) + (c_c * fVcs2);
    qbar = (us_c * fVcd2) + (us * fVud2) + (s   * fVus2) + (s_c * fVcs2);  
  }
  else
  {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }

  LOG("DISSF", pDEBUG) << "Q(x,Q2) = " << q << ", Qbar(x,Q2) = " << qbar;
}
//____________________________________________________________________________
