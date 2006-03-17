//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelNC

\brief    Form Factors for neutrino - free nucleon DIS NC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\ref      E.A.Paschos and J.Y.Yu, Phys.Rev.D 65.033002

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
DISStructureFuncModel("genie::DISStructureFuncModelNC")
{

}
//____________________________________________________________________________
DISStructureFuncModelNC::DISStructureFuncModelNC(string config):
DISStructureFuncModel("genie::DISStructureFuncModelNC", config)
{

}
//____________________________________________________________________________
DISStructureFuncModelNC::~DISStructureFuncModelNC()
{

}
//____________________________________________________________________________
void DISStructureFuncModelNC::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  const Kinematics & kine  = interaction->GetKinematics();
  double x = kine.x();
  if(x<=0. || x>1) {
     LOG("DISSF", pERROR)
                 << "scaling variable x = " << x << "! Can not compute SFs";
     return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the slow rescaling variable and the
  // K factors are applied

  this->CalcPDFs(interaction);

  double u    = fPDF->UpValence() + fPDF->UpSea();
  double ubar = fPDF->UpSea();
  double d    = fPDF->DownValence() + fPDF->DownSea();
  double dbar = fPDF->DownSea();
  double s    = fPDF->Strange();
  double sbar = fPDF->Strange();
  double c    = fPDF->Charm();
  double cbar = fPDF->Charm();

  const InitialState & init_state = interaction->GetInitialState();

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  double F2  = 0.;
  double xF3 = 0.;

  if(isP) {
      F2  = 2*( (kGL2 + kGR2) * (u + c + ubar + cbar) +
                            (kGLprime2 + kGRprime2) * (d + s + dbar + sbar) );
      xF3 = 2*( (kGL2 - kGR2) * (u + c - ubar - cbar) +
                            (kGLprime2 - kGRprime2) * (d + s - dbar - sbar) );

  } else if (isN) {
      F2  = 2*( (kGL2 + kGR2) * (d + c + dbar + cbar) +
                            (kGLprime2 + kGRprime2) * (u + s + ubar + sbar) );
      xF3 = 2*( (kGL2 - kGR2) * (d + c - dbar - cbar) +
                            (kGLprime2 - kGRprime2) * (u + s - ubar - sbar) );
  } else {
     LOG("DISSF", pWARN) << "N type is not handled" << *interaction;
     return;
  }

  // compute nuclear modification factor (for A>1) and the longitudinal
  // structure function (scale-breaking QCD corrections)

  double f  = this->NuclMod(interaction);
  double FL = this->FL(interaction);

  LOG("DISSF", pDEBUG) << "Nucl. Factor = " << f;
  LOG("DISSF", pDEBUG) << "FL = " << FL;

  // compute the structure functions F1-F6
  fF6 = 0.;

  fF3 = f * xF3/x;
  fF2 = f * F2;

  fF1 = 0.5*(fF2-FL)/x; // Callan-Gross holds only if FL=0

  fF5 = fF2/x;          // Albright-Jarlskog relations
  fF4 = 0.;             // Nucl.Phys.B 84, 467 (1975)
}
//____________________________________________________________________________

