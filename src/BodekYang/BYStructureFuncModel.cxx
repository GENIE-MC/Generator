//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModel

\brief    Abstract class.
          Implements part of the DISFormFactorsModelI interface and
          provides some common implementation for concrete Bodek-Yang
          form factor algorithms.

\ref      U.K.Yang and A.Bodek,
          Modeling Deep Inelastic Cross Sections in the Few GeV Region,
          NuINT-01 Proceedings

          U.K.Yang and A.Bodek,
          Parton distributions, d/u and higher twist effect at high x,
          PRL 82, 2467 (1999), hep-ph/9809480,
          PRL 84, 5456 (2000), hep-ph/9912543

          U.K.Yang and A.Bodek,
          Studies of Hugher Twist and Higher Order Effects in NLO and
          NNLO QCD analysis of lepton-nucleon scattering data on F2 and R,
          Eur.Phys.J C13 245, 2000, hep-ex/9908058

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "BodekYang/BYStructureFuncModel.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStructureFuncModel::BYStructureFuncModel() :
DISStructureFuncModel()
{
  this->Init();
}
//____________________________________________________________________________
BYStructureFuncModel::BYStructureFuncModel(const char * param_set):
DISStructureFuncModel(param_set)
{
  this->Init();
}
//____________________________________________________________________________
BYStructureFuncModel::~BYStructureFuncModel()
{

}
//____________________________________________________________________________
void BYStructureFuncModel::Configure(const Registry & config)
{
// Overload Algorithm::Configure() to read the config. registry and set
// private data members.
// DISStructureFuncModel::Configure() creates the owned PDF object that gets
// configured with the specified PDFModelI
// For the ReadBYParams() method see below

  DISStructureFuncModel::Configure(config);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStructureFuncModel::Configure(string param_set)
{
  DISStructureFuncModel::Configure(param_set);
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStructureFuncModel::ReadBYParams(void)
{
// Get the Bodek-Yang model parameters A,B,Csea,Cv1,Cv2 from the config.
// registry and set some private data members so as not to accessing the
// registry at every calculation.
//
  assert( fConfig->Exists("A"  ) );
  assert( fConfig->Exists("B"  ) );
  assert( fConfig->Exists("Cs" ) );
  assert( fConfig->Exists("Cv1") );
  assert( fConfig->Exists("Cv2") );

  fA   = fConfig->GetDouble("A"  );
  fB   = fConfig->GetDouble("B"  );
  fCs  = fConfig->GetDouble("Cs" );
  fCv1 = fConfig->GetDouble("Cv1");
  fCv2 = fConfig->GetDouble("Cv2");

  LOG("BodekYang", pDEBUG) << "Using Bodek-Yang param A   = " << fA;
  LOG("BodekYang", pDEBUG) << "Using Bodek-Yang param B   = " << fB;
  LOG("BodekYang", pDEBUG) << "Using Bodek-Yang param Cs  = " << fCs;
  LOG("BodekYang", pDEBUG) << "Using Bodek-Yang param Cv1 = " << fCv1;
  LOG("BodekYang", pDEBUG) << "Using Bodek-Yang param Cv2 = " << fCv2;
}
//____________________________________________________________________________
void BYStructureFuncModel::Init(void)
{
  fA   = 0;
  fB   = 0;
  fCs  = 0;
  fCv1 = 0;
  fCv2 = 0;
}
//____________________________________________________________________________
bool BYStructureFuncModel::CalculatePDFs(const Interaction * interaction) const
{
  //-- Get the kinematical variables (x,Q2)
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  double x  = sc_params.x();
  double Q2 = this->CalcQ2(interaction);

  //-- Calculate Bodek-Yang modified scaling variable xw
  double xw = x * (Q2 + fB) / (Q2 + fA*x);

  //-- Compute the Bodek-Yang PDFs (GRVLO98 + Bodek-Yang corrections)
  //   The fPDF member was created in DISStructureFuncModel::Configure()
  //    and a BYPDFModel should have been set as the input PDFModelI.
  fPDF->Calculate(xw, Q2);

  //-- scale with sea and valence with the k factors

  double GD   = 1. / TMath::Power(1.+Q2/0.71, 2); // p elastic form factor
  double GD2  = TMath::Power(GD,2);
  double kval = (1.-GD2)*(Q2+fCv2)/(Q2+fCv1);
  double ksea = Q2/(Q2+fCs);

  fPDF->ScaleValence(kval);
  fPDF->ScaleSea(ksea);

  return true;
}
//____________________________________________________________________________
