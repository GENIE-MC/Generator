//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBUtils

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_DB_UTILS_H_
#define _NUVLD_DB_UTILS_H_

#include <cstdio>
#include <iostream>
#include <string>

#include <TGProgressBar.h>
#include <TGTextEdit.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "Utils/StringUtils.h"

using std::string;
using std::cerr;
using std::endl;

using namespace genie::nuvld;

namespace genie {
namespace nuvld {

class DBUtils {

public:

  //___________________________________________________________________________
//  static TSQLResult * send_sql_to_rdbms(const char * sql)
//  {
//    if( is_connected() ) return g_dbase_connection->sql_server()->Query(sql);
//    else  return 0;
//  }
  //___________________________________________________________________________
  static string read_sql_query_from_file(const char * filename)
  {
    // read SQL query
    FILE *fp = fopen(filename, "r");

    if(!fp) {
       cerr << "File " << filename << " could not be read" << endl;
       exit(1);
    }

    char sql[4096] = {};
    fread(sql, 1, 4096, fp);
    fclose(fp);

    string ssql = string(sql);
    return utils::str::FilterString(";", ssql);
  }
  //____________________________________________________________________________
  static void print_sql_result_in_data_viewer(TSQLResult * res)
  {
    //syslog_singleton * syslog = syslog_singleton::instance();
    
    const int nrows = res->GetRowCount();
    //g_data_viewer->AddLine( nuvld_str_utils::concatenate(
    //                                    "number of rows in result: ", nrows) );
    const int nfields = res->GetFieldCount();
    //g_data_viewer->AddLine( nuvld_str_utils::concatenate(
    //                                  "number of fields per row: ", nfields) );

    string * field_name = new string[nfields];

    for (int i = 0; i < nfields; i++)
                                field_name[i] = string( res->GetFieldName(i) );
    TSQLRow * row = 0;

    if(nrows > 0) {

       //double dprogress = 100. / nrows;
       for (int i = 0; i < nrows; i++) {
          row = res->Next();
          //g_data_viewer->AddLine( nuvld_str_utils::concatenate(
          //                               "---------------------- row: ", i) );
          // print all fields
          for (int j = 0; j < nfields; j++) {
              //g_data_viewer->AddLine(nuvld_str_utils::concatenate(
              //             field_name[j].c_str(), " : ", row->GetField(j) ) );
           }
           //syslog->progress_bar()->SetPosition( (i+1)*dprogress );
       }
    }
    //syslog->progress_bar()->SetPosition(0);

    delete [] field_name;      
  }
  //____________________________________________________________________________
    
};

} // nuvld namespace
} // genie namespace

#endif


