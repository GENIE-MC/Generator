//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
 @ Aug 25, 2009 - CA
   Removed redundant versions of ParserUtils.h and ParserStatus.h in favor of
   the ones in $GENIE/Conventions and $GENIE/Utils. Updated code accordingly.

*/
//____________________________________________________________________________ 

#include <cassert>
#include <sstream>
#include <vector>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/SqlQueryBuilder.h"
#include "ValidationTools/NuVld/DBTableFields.h"
#include "ValidationTools/NuVld/DBNuXSecTableFields.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableFields.h"
#include "ValidationTools/NuVld/DBSFTableFields.h"

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
       keys = utils::str::Split(key_list,  ";");

       if(keys.size() == 0) query << " 0 ";
       else {
          query << " (";

          vector<string>::iterator key_iter;

          for(key_iter = keys.begin(); key_iter != keys.end(); ++key_iter) {

             vector<string> key_elem;
             key_elem = utils::str::Split( key_iter->c_str(),  ",");

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
    LOG("NuVld", pINFO) << "Cut List: " << cut_list;

    // add cuts to the SQL query

    if(cut_list.size() > 0) {

        DBTableType_t table_type = query_string.TableType();
        LOG("NuVld", pINFO) << "Table type: " << DBTableType::AsString(table_type);

        vector<string> cuts = utils::str::Split(cut_list,  ";");

        vector<string>::iterator cut_iter;

        switch(table_type) {

        case eDbt_NuXSec:

          // cuts for the CROSS_SECTION table data
          for(cut_iter = cuts.begin(); cut_iter != cuts.end(); ++cut_iter) {

              if(cut_iter->find("Emin") != string::npos) {
                  query << " AND CROSS_SECTION.E > "<< this->CutValue(*cut_iter);
               }
               if(cut_iter->find("Emax") != string::npos) {
                  query << " AND CROSS_SECTION.E < "<< this->CutValue(*cut_iter);
               }
          } // CROSS_SECTION table cuts
          break;

        case eDbt_ElDiffXSec:

          // cuts for the E_DIFF_CROSS_SECTION table data
          for(cut_iter = cuts.begin(); cut_iter != cuts.end(); ++cut_iter) {

              if(cut_iter->find("E_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.E > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("E_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.E < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("EP_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.EP > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("EP_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.EP < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Theta_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Theta > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Theta_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Theta < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Q2_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Q2 > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Q2_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Q2 < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("W2_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.W2 > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("W2_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.W2 < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Nu_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Nu > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Nu_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Nu < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Epsilon_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Epsilon > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Epsilon_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Epsilon < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Gamma_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Gamma > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("Gamma_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.Gamma < "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("x_min") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.x > "
                        << this->CutValue(*cut_iter);
              }
              if(cut_iter->find("x_max") != string::npos) {
                  query << " AND E_DIFF_CROSS_SECTION.x < "
                        << this->CutValue(*cut_iter);
              }
          } // E_DIFF_CROSS_SECTION table cuts
          break;

        case eDbt_SF:

          // cuts for the STRUCTURE_FUNCTION table data
          for(cut_iter = cuts.begin(); cut_iter != cuts.end(); ++cut_iter) {

              if(cut_iter->find("Q2min") != string::npos) {
                  query << " AND STRUCTURE_FUNCTION.Q2 >= "
                        << this->CutValue(*cut_iter);
              } else
              if(cut_iter->find("Q2max") != string::npos) {
                  query << " AND STRUCTURE_FUNCTION.Q2 <= "
                        << this->CutValue(*cut_iter);
              } else
/*
              if(cut_iter->find("x") != string::npos) {
                  query << " AND STRUCTURE_FUNCTION.x > "
                        << this->CutValue(*cut_iter) - 1E-6
                        << " AND STRUCTURE_FUNCTION.x < "
                        << this->CutValue(*cut_iter) + 1E-6;
              } else
*/              
              if(cut_iter->find("xmin") != string::npos) {
                  query << " AND STRUCTURE_FUNCTION.x >= "
                        << this->CutValue(*cut_iter);
              } else
              if(cut_iter->find("xmax") != string::npos) {
                  query << " AND STRUCTURE_FUNCTION.x <= "
                        << this->CutValue(*cut_iter);
              } else

              if(cut_iter->find("R") != string::npos) {

                 query << " AND (";

                 vector<string>::iterator str_iter;
                 vector<string> R = utils::str::Split(
                                         this->CutValueStr(*cut_iter), ",");
                 unsigned int ir=0;
                 for(str_iter = R.begin(); str_iter != R.end(); ++str_iter) {
                    query << " STRUCTURE_FUNCTION.R = \"" 
                          << utils::str::FilterString(" ", *str_iter) << "\"";
                    if(ir++ < R.size()-1) query << " OR ";
                 }
                 query << ")";
              } 
          } // STRUCTURE_FUNCTION table cuts
          break;

        default:
               return "";
        } // db-table-type

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
               fields     = new DBNuXSecTableFields;
               table_name = "CROSS_SECTION";
               break;

      case eDbt_ElDiffXSec:
               fields     = new DBElDiffXSecTableFields;
               table_name = "E_DIFF_CROSS_SECTION";
               break;

      case eDbt_SF:
               fields     = new DBSFTableFields;
               table_name = "STRUCTURE_FUNCTION";
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
               table_name = "CROSS_SECTION";        break;
      case eDbt_ElDiffXSec:
               table_name = "E_DIFF_CROSS_SECTION"; break;
      case eDbt_SF:
               table_name = "STRUCTURE_FUNCTION";   break;
      default:
               return "";
  }
  
  query << ", MEASUREMENT_HEADER WHERE "
        << table_name << ".name = MEASUREMENT_HEADER.name AND "
        << table_name << ".measurement_tag = MEASUREMENT_HEADER.measurement_tag";

  return query.str();        
}
//____________________________________________________________________________
double SqlQueryBuilder::CutValue(string cut_segment)
{
// Extract the cut value from a segment of the cuts sector in DBQueryString
// The cuts sector is a ; separated list of segments CUTS=cut1;cut2;...;cutN
// Each segment is a name=pair string

  vector<string> cut = utils::str::Split(cut_segment, "=");
  assert(cut.size() == 2);
  return atof(cut[1].c_str());
}
//____________________________________________________________________________
string SqlQueryBuilder::CutValueStr(string cut_segment)
{
// Extract the cut value from a segment of the cuts sector in DBQueryString
// The cuts sector is a ; separated list of segments CUTS=cut1;cut2;...;cutN
// Each segment is a name=pair string

  vector<string> cut = utils::str::Split(cut_segment, "=");
  assert(cut.size() == 2);
  return cut[1];
}
//____________________________________________________________________________
