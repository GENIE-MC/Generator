//____________________________________________________________________________
/*!

\namespace  genie::utils::print

\brief      Simple printing utilities

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

\cpright    Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
            All rights reserved.
            For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PRINT_UTILS_H_
#define _PRINT_UTILS_H_

#include <string>

#include <TVector3.h>
#include <TLorentzVector.h>

using std::string;

namespace genie {
namespace utils {

namespace print
{
  string P4AsString      (const TLorentzVector * p);
  string P4AsShortString (const TLorentzVector * p);
  string X4AsString      (const TLorentzVector * x);
  string P3AsString      (const TVector3 * vec);
  string Vec3AsString    (const TVector3 * vec);
  string BoolAsString    (bool b);
  string BoolAsTFString  (bool b);
  string BoolAsIOString  (bool b);
  string BoolAsYNString  (bool b);
  void   PrintBanner     (void);
  string PrintFramedMesg (string mesg);

}      // print namespace
}      // utils namespace
}      // genie namespace

#endif // _PRINT_UTILS_H_
