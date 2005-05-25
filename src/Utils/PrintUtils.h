//____________________________________________________________________________
/*!

\namespace  genie::print_utils

\brief      Simple printing utilities

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#ifndef _PRINT_UTILS_H_
#define _PRINT_UTILS_H_

#include <string>

#include <TVector3.h>
#include <TLorentzVector.h>

using std::string;

namespace genie {

namespace print_utils
{
  string P4AsString   (const TLorentzVector * p);
  string X4AsString   (const TLorentzVector * x);
  string Vec3AsString (const TVector3 * vec);
  string BoolAsString (bool tf);

}      // print_utils namespace

}      // genie namespace

#endif // _PRINT_UTILS_H_
