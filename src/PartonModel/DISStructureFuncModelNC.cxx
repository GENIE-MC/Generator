//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelNC

\brief    Form Factors for neutrino - free nucleon DIS NC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISStructureFuncModelNC.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISStructureFuncModelNC::DISStructureFuncModelNC() :
DISStructureFuncModel()
{
  fName = "genie::DISStructureFuncModelNC";
}
//____________________________________________________________________________
DISStructureFuncModelNC::DISStructureFuncModelNC(const char * param_set):
DISStructureFuncModel(param_set)
{
  fName = "genie::DISStructureFuncModelNC";

  this->FindConfig();
}
//____________________________________________________________________________
DISStructureFuncModelNC::~DISStructureFuncModelNC()
{

}
//____________________________________________________________________________
double DISStructureFuncModelNC::xF1(const Interaction * interaction) const
{
  return DISStructureFuncModel::xF1(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelNC::F2(const Interaction * interaction) const
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

  if(isP)
      F2 = 2*( (kGL2 + kGR2) * (u + c + ubar + cbar) +
                            (kGLprime2 + kGRprime2) * (d + s + dbar + sbar) );

  else if (isN)
      F2 = 2*( (kGL2 + kGR2) * ( d + c + dbar + cbar) +
                            (kGLprime2 + kGRprime2) * (u + s + ubar + sbar) );

  else {
     LOG("DISSF", pWARN) << "N type is not handled" << *interaction;
  }

  return F2;
}
//____________________________________________________________________________
double DISStructureFuncModelNC::xF3(const Interaction * interaction) const
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

  if(isP)
       xF3 = 2* ( (kGL2 - kGR2) * (u + c - ubar - cbar) +
                            (kGLprime2 - kGRprime2) * (d + s - dbar - sbar) );

  else if (isN)
       xF3 = 2* ( (kGL2 - kGR2) * (d + c - dbar - cbar) +
                            (kGLprime2 - kGRprime2) * (u + s - ubar - sbar) );

  else {
     LOG("DISSFs", pWARN) << "N type is not handled" << *interaction;
  }

  return xF3;
}
//____________________________________________________________________________
double DISStructureFuncModelNC::F4(const Interaction * interaction) const
{
  return DISStructureFuncModel::F4(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelNC::xF5(const Interaction * interaction) const
{
  return DISStructureFuncModel::xF5(interaction);
}
//____________________________________________________________________________
double DISStructureFuncModelNC::F6(const Interaction * interaction) const
{
  return DISStructureFuncModel::F6(interaction);
}
//____________________________________________________________________________

