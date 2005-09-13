//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModelCC

\brief    Computes CC vN DIS Form Factors according to the Bodek-Yang model.
          Inherits part of its implemenation from the BYStructureFuncModel
          abstract class.

          Check out BYStructureFuncModel for comments and references.

          Is a concrete implementation of the DISFormFactorsModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#include "BodekYang/BYStructureFuncModelCC.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PDF/PDF.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC() :
BYStructureFuncModel()
{
  fName = "genie::BYStructureFuncModelCC";
}
//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC(const char * param_set):
BYStructureFuncModel(param_set)
{
  fName = "genie::BYStructureFuncModelCC";

  this->FindConfig();
}
//____________________________________________________________________________
BYStructureFuncModelCC::~BYStructureFuncModelCC()
{

}
//____________________________________________________________________________
double BYStructureFuncModelCC::xF1(const Interaction * interaction) const
{
  return BYStructureFuncModel::xF1(interaction);
}
//____________________________________________________________________________
double BYStructureFuncModelCC::F2(const Interaction * interaction) const
{
  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  // Compute F2

  double F2 = 0;

  const InitialState & init_state = interaction->GetInitialState();

  bool isP     = init_state.GetTarget().IsProton();
  bool isN     = init_state.GetTarget().IsNeutron();

  bool isNu    = pdg::IsNeutrino     ( init_state.GetProbePDGCode() );
  bool isNuBar = pdg::IsAntiNeutrino ( init_state.GetProbePDGCode() );

  if      ( isP &&  isNu    )  F2 = 2*(d*kCos8c_2 + s*kSin8c_2 + ubar + cbar);
  else if ( isP &&  isNuBar )  F2 = 2*(u*kCos8c_2 + c*kSin8c_2 + dbar + sbar);
  else if ( isN &&  isNu    )  F2 = 2*(u*kCos8c_2 + s*kSin8c_2 + dbar + cbar);
  else if ( isN &&  isNuBar )  F2 = 2*(d*kCos8c_2 + c*kSin8c_2 + ubar + sbar);
  else {
     LOG("BodekYang", pWARN) << "v/N types are not handled" << *interaction;
  }

  return F2;
}
//____________________________________________________________________________
double BYStructureFuncModelCC::xF3(const Interaction * interaction) const
{
  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  // Compute xF3

  double xF3 = 0;

  const InitialState & init_state = interaction->GetInitialState();

  bool isP     = init_state.GetTarget().IsProton();
  bool isN     = init_state.GetTarget().IsNeutron();

  bool isNu    = pdg::IsNeutrino     ( init_state.GetProbePDGCode() );
  bool isNuBar = pdg::IsAntiNeutrino ( init_state.GetProbePDGCode() );

  if      ( isP && isNu     ) xF3 = 2*(d*kCos8c_2 + s*kSin8c_2 - ubar - cbar);
  else if ( isP &&  isNuBar ) xF3 = 2*(u*kCos8c_2 + c*kSin8c_2 - dbar - sbar);
  else if ( isN &&  isNu    ) xF3 = 2*(u*kCos8c_2 + s*kSin8c_2 - dbar - cbar);
  else if ( isN &&  isNuBar ) xF3 = 2*(d*kCos8c_2 + c*kSin8c_2 - ubar - sbar);
  else {
     LOG("BodekYang", pWARN) << "v/N types are not handled" << *interaction;
  }

  return xF3;
}
//____________________________________________________________________________
double BYStructureFuncModelCC::F4(const Interaction * interaction) const
{
  return BYStructureFuncModel::F4(interaction);
}
//____________________________________________________________________________
double BYStructureFuncModelCC::xF5(const Interaction * interaction) const
{
  return BYStructureFuncModel::xF5(interaction);
}
//____________________________________________________________________________
double BYStructureFuncModelCC::F6(const Interaction * interaction) const
{
  return BYStructureFuncModel::F6(interaction);
}
//____________________________________________________________________________
