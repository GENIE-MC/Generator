//____________________________________________________________________________
/*!

\namespace  genie::string_utils

\brief      Utilities for string manipulation

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    January 12, 2004

*/
//____________________________________________________________________________

#ifndef _STRING_UTILS_H_
#define _STRING_UTILS_H_

#include <string>

using std::string;

namespace genie {

namespace string_utils {

  //! concatenation methods used for NuValidator's TGTextEdit widgets

  const char * Concat(const char * s1, const char * s2,
                                  const char * s3 = 0, const char * s4 = 0);
  const char * Concat(const char * s1, bool b);
  const char * Concat(const char * s1, int n);
  const char * Concat(const char * s1, float x);
  const char * Concat(const char * s1, double x);

} // string_utils namespace
} // genie        namespace

#endif // _STRING_UTILS_H_


