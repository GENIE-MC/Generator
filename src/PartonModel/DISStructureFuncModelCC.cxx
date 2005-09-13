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
void DISStructureFuncModelCC::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the slow rescaling variable and the
  // K factors are applied
  this->CalcPDFs(interaction);

  // Compute q and qbar
  double q    = this -> Q    (interaction);
  double qbar = this -> QBar (interaction);

  if(q<0 || qbar<0) {
     LOG("DISSF", pERROR) << "Negative q and/or q{bar}! Can not compute SFs";
     return;
  }

  const ScatteringParams & sc_params  = interaction->GetScatteringParams();
  double x = sc_params.x();
  if(x<=0.) {
     LOG("DISSF", pERROR)
                 << "scaling variable x = " << x << ". Can not compute SFs";
     return;
  }

  double f = this->NuclMod(interaction);

  fF6 = 0.;
  fF5 = 0.;
  fF4 = 0.;
  fF3 = f * 2*(q-qbar)/x;
  fF2 = f * 2*(q+qbar);
  fF1 = 0.5*fF2/x;
}
//____________________________________________________________________________

