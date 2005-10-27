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
BYStructureFuncModel("genie::BYStructureFuncModelCC")
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC(string config):
BYStructureFuncModel("genie::BYStructureFuncModelCC", config)
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::~BYStructureFuncModelCC()
{

}
//____________________________________________________________________________
void BYStructureFuncModelCC::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the rescaling variable and the
  // K factors are applied
  this->CalcPDFs(interaction);

  // Compute q and qbar
  double q    = this -> Q    (interaction);
  double qbar = this -> QBar (interaction);

  if(q<0 || qbar<0) {
     LOG("BodekYang", pERROR) << "Negative q and/or q{bar}! Can not compute SFs";
     return;
  }

  const Kinematics & kine  = interaction->GetKinematics();
  double x = kine.x();
  if(x<=0.) {
     LOG("BodekYang", pERROR)
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
