//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Provides common implementation for concrete
          DISStructureFuncModelI objects

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModel.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineLimits.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel() :
DISStructureFuncModelI()
{
                     // at each calculation are evaluated at:
  fPDF  = new PDF(); //   x = computed (+/-corrections) scaling var, Q2
  fPDFc = new PDF(); //   x = computed charm slow re-scaling var,    Q2
}
//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel(const char * param_set):
DISStructureFuncModelI(param_set)
{
  fPDF  = new PDF();
  fPDFc = new PDF();
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
  fPDF->SetModel(pdf_model);
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
// Return Q2 from the scattering param objects or, if not set, compute if
// from x,y

  return kine_limits::CalcQ2(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModel::ScalingVar (
                                        const Interaction * interaction) const
{
// This is an abstract class: no model-specific correction
// The scaling variable is set to the normal Bjorken x.
// Override DISStructureFuncModel::ScalingVar() to compute corrections

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  double x  = sc_params.x();

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
  return 1.;
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

  fPDF->Calculate(x, Q2);

  //-- check whether it is above charm threshold
  bool isAbvCh = kine_limits::IsAboveCharmThreshold(interaction, kMc);
  double xc=0.;
  if(isAbvCh) {
    xc = kine_limits::SlowRescalingVar(interaction, kMc);
      fPDF->Calculate(xc, Q2);
  }

  //-- apply the K factors

  double kval = this->KVal(interaction);
  double ksea = this->KSea(interaction);

  fPDF  -> ScaleValence (kval);
  fPDF  -> ScaleSea     (ksea);
  fPDFc -> ScaleValence (kval);
  fPDFc -> ScaleSea     (ksea);
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
  double c    = fPDF  -> Charm();
  double uv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  double us_c = fPDFc -> UpSea();       // ...
  double dv_c = fPDFc -> DownValence(); // ...
  double ds_c = fPDFc -> DownSea();     // ...
  double s_c  = fPDFc -> Strange();     // ...
  double c_c  = fPDFc -> Charm();       // ...

  //-- get CKM constants
  double Vud2 = kVud_2;
  double Vus2 = kVus_2;
  double Vcd2 = kVcd_2;
  double Vcs2 = kVcs_2;

  if (isP && isNu)
  {
     q = (dv+ds)*Vud2 + s*Vus2 + (dv_c+ds_c)*Vcd2 + s_c*Vcs2;
  }
  else if (isP && isNuBar)
  {
     q = (uv+us)*(Vud2+Vus2) + (c_c)*(Vcd2+Vcs2);
  }
  else if (isN && isNu)
  {
     q = (uv+us)*Vud2 + s*Vus2 + (uv_c+us_c)*Vcd2 + s_c*Vcs2;
  }
  else if (isN && isNuBar)
  {
     q = (dv+ds)*(Vud2+Vus2) + (c_c)*(Vcd2+Vcs2);
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
  double uv   = fPDF  -> UpValence();
  double us   = fPDF  -> UpSea();
  double dv   = fPDF  -> DownValence();
  double ds   = fPDF  -> DownSea();
  double s    = fPDF  -> Strange();
  double c    = fPDF  -> Charm();
  double uv_c = fPDFc -> UpValence();   // will be 0 if < charm threshold
  double us_c = fPDFc -> UpSea();       // ...
  double dv_c = fPDFc -> DownValence(); // ...
  double ds_c = fPDFc -> DownSea();     // ...
  double s_c  = fPDFc -> Strange();     // ...
  double c_c  = fPDFc -> Charm();       // ...

  //-- get CKM constants
  double Vud2 = kVud_2;
  double Vus2 = kVus_2;
  double Vcd2 = kVcd_2;
  double Vcs2 = kVcs_2;

  if (isP && isNu)
  {
     qbar = (us)*(Vud2+Vus2) + (c_c)*(Vcd2+Vcs2);
  }
  else if (isP && isNuBar)
  {
     qbar = (ds_c)*(Vcd2) + (ds)*(Vud2) + (s)*(Vus2) + (s_c)*(Vcs2);
  }
  else if (isN && isNu)
  {
     qbar = (ds)*(Vud2+Vus2) + (c_c)*(Vcd2+Vcs2);
  }
  else if (isN && isNuBar)
  {
     qbar = (us_c)*(Vcd2) + (us)*(Vud2) + (s)*(Vus2) + (s_c)*(Vcs2);
  }
  else
  {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }
  return qbar;
}
//____________________________________________________________________________


