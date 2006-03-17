//____________________________________________________________________________
/*!

\class    genie::DISStructureFuncModelCC

\brief    Form Factors for neutrino - free nucleon DIS CC interactions.
          Is a concrete implementation of the DISStructureFuncModelI interface.

\ref      R.Devenish and A.Cooper-Sarkar, Deep Inelastic Scattering, OUP 2004;
          E.A.Paschos and J.Y.Yu, hep-ph/0107261

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
DISStructureFuncModel("genie::DISStructureFuncModelCC")
{

}
//____________________________________________________________________________
DISStructureFuncModelCC::DISStructureFuncModelCC(string config):
DISStructureFuncModel("genie::DISStructureFuncModelCC", config)
{

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

  const Kinematics & kinematics = interaction->GetKinematics();
  double x = kinematics.x();
  if(x<=0. || x>=1) {
     LOG("DISSF", pERROR)
                 << "scaling variable x = " << x << ". Can not compute SFs";
     return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Here all corrections to computing the slow rescaling variable and the
  // K factors are applied
  this->CalcPDFs(interaction);

  // Compute q and qbar
  double q    = this -> Q    (interaction);
  double qbar = this -> QBar (interaction);

  LOG("DISSF", pDEBUG) << "\n" << kinematics;
  LOG("DISSF", pDEBUG) << "q(x,Q2) = " << q << ", q{bar}(x,Q2) = " << qbar;

  if(q<0 || qbar<0) {
     LOG("DISSF", pERROR) << "Negative q and/or q{bar}! Can not compute SFs";
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

  fF3 = f * 2*(q-qbar)/x;
  fF2 = f * 2*(q+qbar);

  fF1 = 0.5*(fF2-FL)/x; // Callan-Gross holds only if FL=0

  fF5 = fF2/x;          // Albright-Jarlskog relations
  fF4 = 0.;             // Nucl.Phys.B 84, 467 (1975)
}
//____________________________________________________________________________

