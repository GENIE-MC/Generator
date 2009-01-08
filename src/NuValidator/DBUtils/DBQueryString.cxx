//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBQueryString

\brief    A collectable string class with encoded information needed to form
          a NuVld database SQL query

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 26, 2005
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <vector>

#include "DBUtils/DBQueryString.h"
#include "XmlParser/ParserUtils.h"

using std::ostringstream;
using std::vector;

using namespace genie::nuvld;

ClassImp(DBQueryString)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  ostream & operator << (ostream & stream, const DBQueryString & query_string)
  {
     stream << query_string.GetString().Data();
     return stream;
  }
 }
}
//____________________________________________________________________________
DBQueryString::DBQueryString()
{

}
//____________________________________________________________________________
DBQueryString::DBQueryString(const DBQueryString & query_string) :
TObjString(query_string.AsString().c_str())
{

}
//____________________________________________________________________________
DBQueryString::DBQueryString(const char * query_string) :
TObjString(query_string)
{

}
//____________________________________________________________________________
DBQueryString::DBQueryString(string query_string) :
TObjString(query_string.c_str())
{

}
//____________________________________________________________________________
DBQueryString::~DBQueryString()
{

}
//____________________________________________________________________________
bool DBQueryString::IsNull(void) const
{
  return false;
}
//____________________________________________________________________________
string DBQueryString::AsString(void) const
{
  string query_str = string( this->GetString().Data() );

  return query_str;
}
//____________________________________________________________________________
string DBQueryString::KeyList(void) const
{
// Extract the list of keys from the DBQueryString
//
  string key_list = this->ExtractFieldValue("KEY-LIST");

  return key_list;  
}
//____________________________________________________________________________
string DBQueryString::CutList(void) const
{
// Extract the list of cuts from the DBQueryString
//
  string cut_list = this->ExtractFieldValue("CUTS");

  return cut_list;
}
//____________________________________________________________________________
DBTableType_t DBQueryString::TableType(void) const
{
  string db_type = this->ExtractFieldValue("DB-TYPE");

  return DBTableType::FromString(db_type);
}
//____________________________________________________________________________
string DBQueryString::ExtractFieldValue(string field_name) const
{
  const char * query_str = this->GetString().Data();

  string field_value = "";

  if( query_str != 0 ) {

    vector<string> selection_elements = ParserUtils::split(query_str,  "$");

    vector<string>::iterator element_iter;

    for(element_iter = selection_elements.begin();
                   element_iter != selection_elements.end(); ++element_iter) {

       if ( element_iter->find(field_name) != string::npos ) {

          vector<string> parts = ParserUtils::split( *element_iter,  ":");
          if(parts.size()==2) field_value = parts[1];
       }
    }// selection string elements
  }// != 0

  return field_value;
}
//____________________________________________________________________________

