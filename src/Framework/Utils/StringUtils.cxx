//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - January 12, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include "Framework/Utils/StringUtils.h"

using std::ostringstream;

//____________________________________________________________________________
string genie::utils::str::TrimSpaces(string input)
{
  if( input.find_first_not_of(" \n") != 0)
       input.erase( 0,  input.find_first_not_of(" \n")  );

  if( input.find_last_not_of(" \n") != input.length() )
       input.erase( input.find_last_not_of(" \n")+1, input.length() );

  return RemoveSuccessiveSpaces(input);
}
//____________________________________________________________________________
string genie::utils::str::IntAsString(int i)
{
  ostringstream os;
  os << i;
  return os.str();
}
//____________________________________________________________________________
vector<string> genie::utils::str::Split(string input, string delimiter)
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
string genie::utils::str::RemoveSuccessiveSpaces(string input)
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
void genie::utils::str::ReplaceStringInPlace(
  string& subject, const string& search, const string& replace)
{
// Searches the string "input" for "search" and replaces this with "replace"
  size_t pos = 0;
  while ((pos = subject.find(search, pos)) != std::string::npos) {
    subject.replace(pos, search.length(), replace);
    // move to end of replaced string
    pos += replace.length();
  }
}
//____________________________________________________________________________
string genie::utils::str::FilterString(string filt_elements, string input)
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
string genie::utils::str::ToUpper(string input)
{
  for(unsigned int i=0; i<input.size(); i++) input[i] = toupper(input[i]);
  return input;
}
//____________________________________________________________________________
string genie::utils::str::ToLower(string input)
{
  for(unsigned int i=0; i<input.size(); i++) input[i] = tolower(input[i]);
  return input;
}
//____________________________________________________________________________

