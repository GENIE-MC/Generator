//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelCC

\brief    Form Factors for neutrino - free nucleon DIS CC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModelCC.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
DISStructureFuncModelCC::DISStructureFuncModelCC() :
DISStructureFuncModel()
{
  fName = "genie::DISStructureFuncModelCC";
}
//____________________________________________________________________________
DISStructureFuncModelCC::DISStructureFuncModelCC(const char * param_set):
DISStructureFuncModel(param_set)
{
  fName = "genie::DISStructureFuncModelCC";

  this->FindConfig();
}
//____________________________________________________________________________
DISStructureFuncModelCC::~DISStructureFuncModelCC()
{

}
//____________________________________________________________________________
double DISStructureFuncModelCC::xF1(const Interaction * interaction) const
{
  return DISStructureFuncModel::xF1(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelCC::F2(const Interaction * interaction) const
{
  //this->CalcPDFs(interaction);

  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  double F2 = 0;

  const InitialState & init_state = interaction->GetInitialState();

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  bool isNu    = pdg::IsNeutrino     ( init_state.GetProbePDGCode() );
  bool isNuBar = pdg::IsAntiNeutrino ( init_state.GetProbePDGCode() );

  if      ( isP && isNu     )  F2  = 2*(d+s+ubar+cbar);
  else if ( isP &&  isNuBar )  F2  = 2*(u+c+dbar+sbar);
  else if ( isN &&  isNu    )  F2  = 2*(u+s+dbar+cbar);
  else if ( isN &&  isNuBar )  F2  = 2*(d+c+ubar+sbar);
  else {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }

  return F2;
}
//____________________________________________________________________________
double DISStructureFuncModelCC::xF3(const Interaction * interaction) const
{
  //this->CalcPDFs(interaction);

  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  double xF3 = 0;

  const InitialState & init_state = interaction->GetInitialState();

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  bool isNu    = pdg::IsNeutrino     ( init_state.GetProbePDGCode() );
  bool isNuBar = pdg::IsAntiNeutrino ( init_state.GetProbePDGCode() );

  if      ( isP && isNu     )  xF3  = 2*(d+s-ubar-cbar);
  else if ( isP &&  isNuBar )  xF3  = 2*(u+c-dbar-sbar);
  else if ( isN &&  isNu    )  xF3  = 2*(u+s-dbar-cbar);
  else if ( isN &&  isNuBar )  xF3  = 2*(d+c-ubar-sbar);
  else {
     LOG("DISSF", pWARN) << "v/N types are not handled" << *interaction;
  }

  return xF3;
}
//____________________________________________________________________________
double DISStructureFuncModelCC::F4(const Interaction * interaction) const
{
  return DISStructureFuncModel::F4(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelCC::xF5(const Interaction * interaction) const
{
  return DISStructureFuncModel::xF5(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelCC::F6(const Interaction * interaction) const
{
  return DISStructureFuncModel::F6(interaction);
}
//____________________________________________________________________________


