//_____________________________________________________________________________
/*!

\class    genie::nuvld::ParserUtils

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#ifndef _PARSER_UTILS_H_
#define _PARSER_UTILS_H_

#include <sstream>
#include <string>
#include <iostream>
#include <vector>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

using std::ostringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

namespace genie {
namespace nuvld {

class ParserUtils {

public:

 //__________________________________________________________________________
 static string ParserUtils::trim_spaces(xmlChar * xml_string)
 {
   // trim the leading/trailing empty spaces from the of the parsed xml string
   // like in
   //
   // "      I am a string with lots of spaces      " ---->
   //                                  "I am a string with lots of spaces"
   //
   // In this method, "\n" is treated as 'empty space' so as to trim not only empty
   // spaces in the line that contains the string but also all leading and trailing 
   // empty lines

   string std_string  = string( (const char *) xml_string );

   if( std_string.find_first_not_of(" \n") != 0) 
        std_string.erase( 0,  std_string.find_first_not_of(" \n")-1  );
	
   if( std_string.find_last_not_of(" \n") != std_string.length() ) 	
        std_string.erase( std_string.find_last_not_of(" \n")+1, std_string.length() ); 

   return remove_successive_spaces(std_string);
 }
 //__________________________________________________________________________
 static string ParserUtils::int_as_string(int i)
 {
   ostringstream os;   os << i;    return os.str();
 }
 //__________________________________________________________________________
 static vector<string> ParserUtils::split(string str, string delimiter)
 {
  // split a string of 'delimiter' separated values and return each string 
  // part as a vector<string> element.

   vector<string> string_parts;
  
   while(str.find_first_of(delimiter) < str.length()) {

     string substring = str.substr(0, str.find_first_of(delimiter));
     
     if(substring.size() > 0) string_parts.push_back( substring );

     str = str.substr(str.find_first_of(delimiter)+1, str.length());
   }

   if(str.size() > 0)  string_parts.push_back(str);
   return string_parts;
 }
 //__________________________________________________________________________
 static string ParserUtils::remove_successive_spaces(string str)
 {
  // this method trims white space that may be found within an expression...
  // eg. "stop    pressing     the spacebar" -> "stop pressing the spacebar"

   string trimmed = str;

   string::size_type pos = 0;

   while( (pos = trimmed.find_first_of(" ", pos)) != string::npos ) {
      if( trimmed[++pos] == ' ' ) trimmed.erase(pos--,1);
      //cout << trimmed << endl;
   }

   return filter_string("\n",trimmed);
 }
 //__________________________________________________________________________
 static string ParserUtils::get_attribute(
                                        xmlNodePtr xml_cur, string attr_name)
 {
   xmlChar *  xml_string = xmlGetProp( 
                             xml_cur, (const xmlChar *) attr_name.c_str() );
   string     std_string = ParserUtils::trim_spaces(xml_string);
  
   xmlFree(xml_string);
  
   return std_string;
 }
 //__________________________________________________________________________
 static string ParserUtils::filter_string(string filt_elements, string input)
 {
  // filter out from 'input' all characters that can be found in 'filt_elements'

   string filtered_string = input;
   string::size_type pos = 0;

   while( (pos = filtered_string.find_first_of(
                       filt_elements, pos)) != string::npos )
                                                filtered_string.erase(pos,1);

   return filtered_string;
 }
 //__________________________________________________________________________

};

} // nuvld namespace
} // genie namespace

#endif // _PARSER_UTILS_H_

