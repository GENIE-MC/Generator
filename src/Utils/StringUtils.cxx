//____________________________________________________________________________
/*!

\namespace  genie::string_utils

\brief      Utilities for string manipulation

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    January 12, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include "Utils/StringUtils.h"

using std::ostringstream;

//____________________________________________________________________________
const char * genie::string_utils::Concat(
           const char * s1, const char * s2, const char * s3, const char * s4)
{
  ostringstream sstr;
  sstr << s1 << s2;
  
  if(s3) sstr << s3;
  if(s4) sstr << s4;
    
  return sstr.str().c_str();
}
//____________________________________________________________________________
const char * genie::string_utils::Concat(const char * s1, bool b)
{
    ostringstream sstr;
    sstr << s1;

    if(b) sstr << "TRUE";
    else  sstr << "FALSE";

    return sstr.str().c_str();
}
//____________________________________________________________________________
const char * genie::string_utils::Concat(const char * s1, int n)
{
  ostringstream sstr;
  sstr << s1 << n;
  return sstr.str().c_str();
}
//____________________________________________________________________________
const char * genie::string_utils::Concat(const char * s1, double x)
{
  ostringstream sstr;
  sstr << s1 << x;
  return sstr.str().c_str();
}    
//____________________________________________________________________________
const char * genie::string_utils::Concat(const char * s1, float x)
{
  ostringstream sstr;
  sstr << s1 << x;
  return sstr.str().c_str();
}    
//____________________________________________________________________________


