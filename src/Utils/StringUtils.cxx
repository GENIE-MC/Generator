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
string genie::string_utils::TrimSpaces(string input)
{
  if( input.find_first_not_of(" \n") != 0)
       input.erase( 0,  input.find_first_not_of(" \n")  );

  if( input.find_last_not_of(" \n") != input.length() )
       input.erase( input.find_last_not_of(" \n")+1, input.length() );

  return RemoveSuccessiveSpaces(input);
}
//____________________________________________________________________________
string genie::string_utils::IntAsString(int i)
{
  ostringstream os;
  os << i;
  return os.str();
}
//____________________________________________________________________________
vector<string> genie::string_utils::Split(string input, string delimiter)
{
// split a string of 'delimiter' separated values and return each string
// part as a vector<string> element.

  vector<string> string_parts;

  while(input.find_first_of(delimiter) < input.length()) {

    string_parts.push_back( input.substr(0, input.find_first_of(delimiter)) );
    input = input.substr(input.find_first_of(delimiter)+1, input.length());
  }
  string_parts.push_back(input);
  return string_parts;
}
//____________________________________________________________________________
string genie::string_utils::RemoveSuccessiveSpaces(string input)
{
// this method trims white space that may be found within an expression...
// eg. "stop    pressing     the spacebar" -> "stop pressing the spacebar"

  string trimmed = input;

  string::size_type pos = 0;

  while( (pos = trimmed.find_first_of(" ", pos)) != string::npos ) {
     if( trimmed[++pos] == ' ' ) trimmed.erase(pos--,1);
  }
  return FilterString("\n",trimmed);
}
//____________________________________________________________________________
string genie::string_utils::FilterString(string filt_elements, string input)
{
// filter out from 'input' all characters that can be found in 'filt_elements'

  string filtered_string = input;
  string::size_type pos = 0;

  while( (pos = filtered_string.find_first_of(
                       filt_elements, pos)) != string::npos )
                                                filtered_string.erase(pos,1);
  return filtered_string;
}
//____________________________________________________________________________
string genie::string_utils::ToUpper(string input)
{
  for(unsigned int i=0; i<input.size(); i++) input[i] = toupper(input[i]);
  return input;
}
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


