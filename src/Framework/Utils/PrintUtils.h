//____________________________________________________________________________
/*!

\namespace  genie::utils::print

\brief      Simple printing utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 06, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
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
  void   PrintBanner     (string filename, UInt_t wait_msec);
  string PrintFramedMesg (string mesg, unsigned int nl=1, const char f='*');

}      // print namespace
}      // utils namespace
}      // genie namespace

#endif // _PRINT_UTILS_H_
