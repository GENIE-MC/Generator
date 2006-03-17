//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Provides common implementation for concrete
          DISStructureFuncModelI objects

\ref      For a discussion of DIS SF see for eaxample E.A.Paschos and J.Y.Yu, 
          Phys.Rev.D 65.033002 and R.Devenish and A.Cooper-Sarkar, OUP 2004.
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModel.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"

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
  this->ConfigPDF();
}
//____________________________________________________________________________
void DISStructureFuncModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->ConfigPDF();
}
//____________________________________________________________________________
void DISStructureFuncModel::ConfigPDF(void)
{
  const PDFModelI * pdf_model =
                  dynamic_cast<const PDFModelI *>
                             (this->SubAlg("pdf-alg-name", "pdf-param-set"));
  fPDF  -> SetModel(pdf_model);
  fPDFc -> SetModel(pdf_model);
}
//____________________________________________________________________________
void DISStructureFuncModel::InitPDF(void)
{
                     // at each calculation are evaluated at:
  fPDF  = new PDF(); //   x = computed (+/-corrections) scaling var, Q2
  fPDFc = new PDF(); //   x = computed charm slow re-scaling var,    Q2

  fQ2min = 0; // concrete DISStructureFuncModelI can set it at their config
}
//____________________________________________________________________________
void DISStructureFuncModel::Calculate(const Interaction * ) const
{
  // this is an abstract class - just reset them
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;
}
//____________________________________________________________________________
double DISStructureFuncModel::Q2(const Interaction * interaction) const
{
// Return Q2 from the kinematics or, if not set, compute it from x,y

  return utils::kinematics::CalcQ2(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModel::ScalingVar (
                                        const Interaction * interaction) const
{
// This is an abstract class: no model-specific correction
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  const Kinematics & kine = interaction->GetKinematics();
  double x = kine.x();
  return x;
}
//____________________________________________________________________________
double DISStructureFuncModel::KSea (const Interaction * ) const
{
// This is an abstract class: no model-specific correction
// The sea PDF scaling variable is set to 1
// Override DISStructureFuncModel::KSea() to compute corrections

  return 1.;
}
//____________________________________________________________________________
double DISStructureFuncModel::KVal (const Interaction * ) const
{
// This is an abstract class: no model-specific correction
// The valence PDF scaling variable is set to 1
// Override DISStructureFuncModel::KVal() to compute corrections

  return 1.;
}
//____________________________________________________________________________
double DISStructureFuncModel::NuclMod(const Interaction * ) const
{
// Nuclear modification to Fi

  return 1.;
}
//____________________________________________________________________________
double DISStructureFuncModel::FL(const Interaction * ) const
{
// Longitudinal structure function - If set violates the Callan-Gross relation

  return 0.;
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
  bool isAbvCh = utils::kinematics::IsAboveCharmThreshold(interaction, kMc);
  if(isAbvCh) {
    // compute PDFs at (xi,Q2)
    double xc = utils::kinematics::SlowRescalingVar(interaction, kMc);    
    if(xc>0 && xc<1) fPDFc->Calculate(xc, Q2);
  }

  //-- apply the K factors
  double kval = this->KVal(interaction);
  double ksea = this->KSea(interaction);

  fPDF->ScaleValence (kval);
  fPDF->ScaleSea     (ksea);
  if(isAbvCh) {
    fPDFc->ScaleValence (kval);
    fPDFc->ScaleSea     (ksea);
  }
}
//____________________________________________________________________________
double DISStructureFuncModel::Q(const Interaction * interaction) const
{
  double q  = -1;

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
  //double c    = fPDF  -> Charm();
  double uv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  double us_c = fPDFc -> UpSea();       // ...
  double dv_c = fPDFc -> DownValence(); // ...
  double ds_c = fPDFc -> DownSea();     // ...
  double s_c  = fPDFc -> Strange();     // ...
  double c_c  = fPDFc -> Charm();       // ...

  // Rules of thumb for computing Q and QBar
  // ---------------------------------------
  // - For W+ exchange use: -1/3|e| quarks and -2/3|e| antiquarks
  // - For W- exchange use:  2/3|e| quarks and  1/3|e| antiquarks
  // - For each qi -> qj transition multiply with the (ij CKM element)^2
  // - Use isospin symmetry to get neutro's u,d from proton's u,d
  //    -- neutron d = proton u
  //    -- neutron u = proton d
  // - Use u = usea + uvalence. Same for d
  // - For s,c use q=qbar
  // - For t,b use q=qbar=0

  if (isP && isNu)
  {
     q = (dv+ds)*kVud_2 + s*kVus_2 + (dv_c+ds_c)*kVcd_2 + s_c*kVcs_2;
  }
  else if (isP && isNuBar)
  {
     q = (uv+us)*(kVud_2+kVus_2) + (c_c)*(kVcd_2+kVcs_2);
  }
  else if (isN && isNu)
  {
     q = (uv+us)*kVud_2 + s*kVus_2 + (uv_c+us_c)*kVcd_2 + s_c*kVcs_2;
  }
  else if (isN && isNuBar)
  {
     q = (dv+ds)*(kVud_2+kVus_2) + (c_c)*(kVcd_2+kVcs_2);
  }
  else
  {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }
  return q;
}
//____________________________________________________________________________
double DISStructureFuncModel::QBar(const Interaction * interaction) const
{
  double qbar  = -1;

  const InitialState & init_state = interaction->GetInitialState();

  int nuc_pdgc = init_state.GetTarget().StruckNucleonPDGCode();
  int nu_pdgc  = init_state.GetProbePDGCode();
  bool isP     = pdg::IsProton       ( nuc_pdgc );
  bool isN     = pdg::IsNeutron      ( nuc_pdgc );
  bool isNu    = pdg::IsNeutrino     ( nu_pdgc  );
  bool isNuBar = pdg::IsAntiNeutrino ( nu_pdgc  );

  //-- get PDFs [should have been computed by calling in CalcPDFs() first]
  //double uv   = fPDF  -> UpValence();
  double us   = fPDF  -> UpSea();
  //double dv   = fPDF  -> DownValence();
  double ds   = fPDF  -> DownSea();
  double s    = fPDF  -> Strange();
  //double c    = fPDF  -> Charm();
  //double uv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  double us_c = fPDFc -> UpSea();       // ...
  //double dv_c = fPDFc -> DownValence(); // ...
  double ds_c = fPDFc -> DownSea();     // ...
  double s_c  = fPDFc -> Strange();     // ...
  double c_c  = fPDFc -> Charm();       // ...

  if (isP && isNu)
  {
     qbar = (us)*(kVud_2+kVus_2) + (c_c)*(kVcd_2+kVcs_2);
  }
  else if (isP && isNuBar)
  {
     qbar = (ds_c)*(kVcd_2) + (ds)*(kVud_2) + (s)*(kVus_2) + (s_c)*(kVcs_2);
  }
  else if (isN && isNu)
  {
     qbar = (ds)*(kVud_2+kVus_2) + (c_c)*(kVcd_2+kVcs_2);
  }
  else if (isN && isNuBar)
  {
     qbar = (us_c)*(kVcd_2) + (us)*(kVud_2) + (s)*(kVus_2) + (s_c)*(kVcs_2);
  }
  else
  {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }
  return qbar;
}
//____________________________________________________________________________


