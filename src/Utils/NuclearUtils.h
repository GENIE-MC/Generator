//____________________________________________________________________________
/*!

\namespace  genie::utils::nuclear

\brief      Simple nuclear physics empirical formulas

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_UTILS_H_
#define _NUCLEAR_UTILS_H_

namespace genie {

class Target;

namespace utils {

namespace nuclear
{
  double BindEnergy            (const Target & target);
  double BindEnergyPerNucleon  (const Target & target);
  double BindEnergyLastNucleon (const Target & target);
  double Radius                (const Target & target);

} // nuclear namespace
} // utils   namespace
} // genie   namespace

#endif // _NUCL_UTILS_H_
