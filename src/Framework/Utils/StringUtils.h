//____________________________________________________________________________
/*!

\namespace  genie::utils::str

\brief      Utilities for string manipulation

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    January 12, 2004

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org            
*/
//____________________________________________________________________________

#ifndef _STRING_UTILS_H_
#define _STRING_UTILS_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace genie {
namespace utils {

namespace str
{
  string         TrimSpaces             (string input);
  string         IntAsString            (int i);
  vector<string> Split                  (string input, string delim);
  string         RemoveSuccessiveSpaces (string input);
  void           ReplaceStringInPlace   (string& subject, const string& search, const string& replace);
  string         FilterString           (string filt, string input);
  string         ToUpper                (string input);
  string         ToLower                (string input);

  template<class T>
    bool Convert( const vector<std::string> & input, std::vector<T> & v ) ;


} // str    namespace
} // utils  namespace
} // genie  namespace

#ifndef __CINT__  // don't even try for ROOT 5
#include "Framework/Utils/StringUtils.icc"
#endif

#endif // _STRING_UTILS_H_
