//____________________________________________________________________________
/*!

\namespace  genie::utils::nuclear

\brief      Simple nuclear physics empirical formulas

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Interaction/Target.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::nuclear::BindEnergy(const Target & target)
{
// Compute the average binding energy (in GeV) using the semi-empirical
// formula from Wapstra (Handbuch der Physik, XXXVIII/1)

  if(!target.IsNucleus()) return 0;

  double a = 15.835;
  double b =  18.33;
  double s =  23.20;
  double d =   0.714;

  double delta = 0;                       /*E-O*/
  if (target.IsOddOdd()  ) delta =  11.2; /*O-O*/
  if (target.IsEvenEven()) delta = -11.2; /*E-E*/

  double N = (double) target.N();
  double Z = (double) target.Z();
  double A = (double) target.A();

  double BE =  a * A
             - b * TMath::Power(A,0.667)
             - s * TMath::Power(N-Z,2.0)/A
             - d * TMath::Power(Z,2.0)/TMath::Power(A,0.333)
             - delta / TMath::Sqrt(A); // MeV

  return (1e-3 * BE); // GeV
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergyPerNucleon(const Target & target)
{
// Compute the average binding energy per nucleon (in GeV)

  if(!target.IsNucleus()) return 0;

  return (utils::nuclear::BindEnergy(target) / target.A());
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergyLastNucleon(const Target & target)
{
// Compute the binding for the most loose nucleon (in GeV)

  if(!target.IsNucleus()) return 0;

  //-- temporarily, return the binding energy per nucleon rather than the
  //   separation energy of the last nucleon
   return (utils::nuclear::BindEnergy(target) / target.A());
}
//___________________________________________________________________________
double genie::utils::nuclear::Radius(const Target & target)
{
// Compute the nuclear radius (in GeV^-1)
// (to get it in F, for example, use: Radius(target)/genie::units::fermi)

  if(!target.IsNucleus()) return 0;

  double A  = (double) target.A();
  double Rn = kNucRo * TMath::Power(A, 0.3333); // in GeV^-1

  return Rn;
}
//___________________________________________________________________________

