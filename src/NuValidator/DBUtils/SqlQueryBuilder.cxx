//_____________________________________________________________________________
/*!

\class    genie::nuvld::SqlQueryBuilder

\brief    Utility class used by the DBI to load DBTable<T>'s from the RDBMS

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <vector>

#include "DBUtils/SqlQueryBuilder.h"
#include "DBUtils/DBTableFields.h"
#include "DBUtils/vXSecTableFields.h"
#include "DBUtils/eDiffXSecTableFields.h"
#include "Messenger/Messenger.h"
#include "XmlParser/ParserUtils.h"

using std::ostringstream;
using std::vector;

using namespace genie::nuvld;

//____________________________________________________________________________
SqlQueryBuilder::SqlQueryBuilder()
{

}
//____________________________________________________________________________
SqlQueryBuilder::~SqlQueryBuilder()
{

}
//____________________________________________________________________________
string SqlQueryBuilder::FormQuery(const DBQueryString & query_string)
{
  ostringstream query;
  
  // add the field names for the table to be queried

  query << this->AddTableFields(query_string);
  
  // add a user-conditional (keys/cuts) to the sql query & perform the join

  ostringstream conditional;

  conditional << this->AddKeyList(query_string) << this->AddCutList(query_string);

  bool have_conditional = conditional.str().length() > 0;
  
  if(have_conditional) query << this->MakeJoin(query_string);

  if(have_conditional) query << " AND " << conditional.str() << ";";
  else                 query << ";";


  SLOG("NuVld", pDEBUG) << "SQL Query = " << query.str();

  return query.str();
}
//____________________________________________________________________________
string SqlQueryBuilder::AddKeyList(const DBQueryString & query_string)
{
  unsigned int i = 0;

  ostringstream query;

  if( ! query_string.IsNull() ) {

    string key_list = query_string.KeyList();

    // add key-list to the SQL query

    if( strcmp(key_list.c_str(),"*") != 0 )  {

       // the key list is formated as: name1,tag1;name2,tag2;name3,tag3;...
       // where each (name,tag) pair is unique dbase key

       vector<string> keys;
       keys = ParserUtils::split(key_list,  ";");

       if(keys.size() == 0) query << " 0 ";
       else {
          query << " (";

          vector<string>::iterator key_iter;

          for(key_iter = keys.begin(); key_iter != keys.end(); ++key_iter) {

             vector<string> key_elem;
             key_elem = ParserUtils::split( key_iter->c_str(),  ",");

             assert(key_elem.size() == 2);

             string name = key_elem[0];
             string mtag = key_elem[1];

             query
                   << "( MEASUREMENT_HEADER.name = \"" << name << "\" AND "
                   << "  MEASUREMENT_HEADER.measurement_tag = \"" << mtag << "\" )";

             // build composite boolean expression - do not add 'OR' after the last key
             if(i++ < keys.size()-1 ) query << " OR ";

          } // number of keys

          query << ") ";
       }

    } else query << " 1 ";

  } else query << " AND 0 ";

  return query.str();
}
//____________________________________________________________________________
string SqlQueryBuilder::AddCutList(const DBQueryString & query_string)
{
  ostringstream query;

  if( ! query_string.IsNull() ) {

    string cut_list = query_string.CutList();

    // add cuts to the SQL query

    if(cut_list.size() > 0) {

       vector<string> cuts = ParserUtils::split(cut_list,  ";");

       vector<string>::iterator cut_iter;

       for(cut_iter = cuts.begin(); cut_iter != cuts.end(); ++cut_iter) {

          if(cut_iter->find("Emin") != string::npos) {

               vector<string> cut_parts = ParserUtils::split(*cut_iter, "=");
               assert(cut_parts.size() == 2);

               query << " AND CROSS_SECTION.E > "
                                              << atof( cut_parts[1].c_str() );
          }

          if(cut_iter->find("Emax") != string::npos) {

               vector<string> cut_parts = ParserUtils::split(*cut_iter, "=");
               assert(cut_parts.size() == 2);

               query << " AND CROSS_SECTION.E < "
                                              << atof( cut_parts[1].c_str() );
          }

       } // cuts
    } // >0 cut-list

  } else query << " AND 0 ";

  return query.str();
}
//____________________________________________________________________________
string SqlQueryBuilder::AddTableFields(const DBQueryString & query_string)
{
  ostringstream query;

  DBTableType_t table_type = query_string.TableType();

  SLOG("NuVld", pINFO) << "Table type: " << DBTableType::AsString(table_type);

  string table_name;
  const DBTableFields * fields = 0;

  switch(table_type) {

      case eDbt_NuXSec:
               fields     = new vXSecTableFields;
               table_name = "CROSS_SECTION";
               break;

      case eDbt_ElDiffXSec:
               fields     = new eDiffXSecTableFields;
               table_name = "E_DIFF_CROSS_SECTION";
               break;
      default:
               return "";
  }

  query << "SELECT ";
  
  for(vector<string>::size_type i = 0; i < fields->NFields(); i++) {
       query << table_name << "." << fields->FieldName(i);
       if(i < fields->NFields() - 1) query << ", ";
  }
  
  query << " FROM " << table_name << " ";

  delete fields;
  
  return query.str();
}
//____________________________________________________________________________
string SqlQueryBuilder::MakeJoin(const DBQueryString & query_string)
{
  ostringstream query;

  DBTableType_t table_type = query_string.TableType();

  string table_name;

  switch(table_type) {
      case eDbt_NuXSec:
               table_name = "CROSS_SECTION"; break;
      case eDbt_ElDiffXSec:
               table_name = "E_DIFF_CROSS_SECTION"; break;
      default:
               return "";
  }
  
  query << ", MEASUREMENT_HEADER WHERE "
        << table_name << ".name = MEASUREMENT_HEADER.name AND "
        << table_name << ".measurement_tag = MEASUREMENT_HEADER.measurement_tag";

  return query.str();        
}
//____________________________________________________________________________

