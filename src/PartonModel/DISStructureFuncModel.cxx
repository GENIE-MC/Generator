//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModel

\brief    Abstract base class. Implements the DISFormFactorsModelI interface
          but can not be instantiated. Its mere role of existence is to factor
          out common implementation from concrete implementations like
          the PartonModelNC, PartonModelCCAboveCharmThr and
          PartonModelCCBelowCharmThr.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModel.h"
#include "PDF/PDFModelI.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel() :
DISStructureFuncModelI()
{
  fPDF = new PDF();
}
//____________________________________________________________________________
DISStructureFuncModel::DISStructureFuncModel(const char * param_set):
DISStructureFuncModelI(param_set)
{
  fPDF = new PDF();
}
//____________________________________________________________________________
DISStructureFuncModel::~DISStructureFuncModel()
{
  delete fPDF;
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
//void DISStructureFuncModel::CalcPDFs(const Interaction * interaction) const
bool DISStructureFuncModel::CalculatePDFs(const Interaction * interaction) const
{
  //-- get scattering parameters
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  double x  = sc_params.x();
  double Q2 = this->CalcQ2(interaction);

  //-- ask PDF obj. to delegate request for calculation to its attached model
  fPDF->Calculate(x, Q2);

  return true;
}
//____________________________________________________________________________
/*
void DISStructureFuncModel::CalcPDFs(double x, double Q2) const
{
  LOG("PartonModel", pDEBUG) << "Computing PDFs";

  //-- ask PDF obj. to delegate request for calculation to its attached model
  fPDF->Calculate(x, Q2);

  LOG("PartonModel", pDEBUG) << *fPDF;
}*/
//____________________________________________________________________________
double DISStructureFuncModel::CalcQ2(const Interaction * interaction) const
{
  const ScatteringParams & sc_params = interaction -> GetScatteringParams();
  const InitialState &    init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double Q2 = 0;

  if( sc_params.Exists("Q2") ) Q2 = sc_params.Q2();
  else if ( sc_params.Exists("y") ) {

    double x  = sc_params.x();
    double y  = sc_params.y();
    double Ev = p4->Energy();
    double Mn = init_state.GetTarget().StruckNucleonMass();

    Q2 = 2*Mn*Ev*x*y;
  }
  delete p4;

  LOG("PartonModel", pDEBUG) << "Q^2 = " << Q2;

  return Q2;
}
//____________________________________________________________________________





