//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModel

\brief    Abstract class. Provides common implementation for concrete
          DISStructureFuncModelI objects computing the Bodek Yang structure
          functions.

\ref      hep-ph/0411202

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
BYStructureFuncModel::BYStructureFuncModel(string name) :
DISStructureFuncModel(name)
{

}
//____________________________________________________________________________
BYStructureFuncModel::BYStructureFuncModel(string name, string config):
DISStructureFuncModel(name, config)
{

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

  this->Init();
  this->ReadBYParams();
}
//____________________________________________________________________________
void BYStructureFuncModel::Configure(string param_set)
{
  DISStructureFuncModel::Configure(param_set);

  this->Init();
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
double BYStructureFuncModel::ScalingVar(const Interaction * interaction) const
{
// Overrides DISStructureFuncModel::ScalingVar() to compute the BY scaling var

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  double x  = sc_params.x();
  double Q2 = this->Q2(interaction);

  double xw = x * (Q2 + fB) / (Q2 + fA*x);

  return xw;
}
//____________________________________________________________________________
double BYStructureFuncModel::KSea(const Interaction * interaction) const
{
// Overrides DISStructureFuncModel::KSea() to compute the BY KSea factor

  double Q2   = this->Q2(interaction);

  double ksea = Q2/(Q2+fCs);
  return ksea;
}
//____________________________________________________________________________
double BYStructureFuncModel::KVal(const Interaction * interaction) const
{
// Overrides DISStructureFuncModel::KSea() to compute the BY KVal factor

  double Q2   = this->Q2(interaction);
  double GD   = 1. / TMath::Power(1.+Q2/0.71, 2); // p elastic form factor
  double GD2  = TMath::Power(GD,2);

  double kval = (1.-GD2)*(Q2+fCv2)/(Q2+fCv1);
  return kval;
}
//____________________________________________________________________________
