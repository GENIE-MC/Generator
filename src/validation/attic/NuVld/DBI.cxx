//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 25, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TSQLResult.h>
#include <TSQLRow.h>

#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBXmlUploader.h"
#include "ValidationTools/NuVld/SqlQueryBuilder.h"
#include "Messenger/Messenger.h"

using namespace genie::nuvld;

//______________________________________________________________________________
DBI::DBI()
{
  fSqlServer = 0;
}
//______________________________________________________________________________
DBI::DBI(TSQLServer * sql_server)
{
  fSqlServer = sql_server;
}
//______________________________________________________________________________
DBI::~DBI()
{

}
//______________________________________________________________________________
DBStatus_t DBI::UploadXML(const XmlDataSet & data)
{
  if( this->HaveConnection() ) {
    
    DBXmlUploader xml_uploader(fSqlServer);
  
    DBStatus_t status = xml_uploader.Upload(data);

    return status;    
  } 

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
DBStatus_t DBI::FillTable(
  DBTable<DBNuXSecTableRow> * table, const DBQueryString & query_string)
{
  if( this->HaveConnection() ) {

     SLOG("NuVld", pINFO)
       << "Forming SQL-Query for DBQueryString: " << query_string;
  
     SqlQueryBuilder qbld;
      
     string sql_query = qbld.FormQuery(query_string);
     SLOG("NuVld", pINFO) << "Formed SQL-Query: " << sql_query;

     TSQLResult * res = fSqlServer->Query( sql_query.c_str() );

     const int nrows = res->GetRowCount();

     SLOG("NuVld", pINFO) << "The Query returned: " << nrows << " rows";

     if(nrows == 0) return eDbu_EMPTY_TABLE;
     
     TSQLRow *       row      = 0;
     DBNuXSecTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new DBNuXSecTableRow(row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;
  
     this  -> AddDBMeasurementIdList (table);
     table -> SetQueryString       (query_string);

     return eDbu_OK;
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
DBStatus_t DBI::FillTable(
      DBTable<DBElDiffXSecTableRow> * table, const DBQueryString & query_string)
{
  if( this->HaveConnection() ) {

     SLOG("NuVld", pINFO)
                     << "Forming SQL-Query for DBQueryString: " << query_string;

     SqlQueryBuilder qbld;

     string sql_query = qbld.FormQuery(query_string);

     SLOG("NuVld", pINFO) << "Formed SQL-Query: " << sql_query;

     TSQLResult * res = fSqlServer->Query( sql_query.c_str() );

     const int nrows = res->GetRowCount();
     SLOG("NuVld", pINFO) << "The Query returned: " << nrows << " rows";

     if(nrows == 0) return eDbu_EMPTY_TABLE;
     
     TSQLRow * row = 0;
     DBElDiffXSecTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new DBElDiffXSecTableRow (row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;

     this  -> AddDBMeasurementIdList (table);
     table -> SetQueryString       (query_string);
     
     return eDbu_OK;     
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
DBStatus_t DBI::FillTable(
   DBTable<DBSFTableRow> * table, const DBQueryString & query_string)
{
  if( this->HaveConnection() ) {

     SLOG("NuVld", pINFO)
       << "Forming SQL-Query for DBQueryString: " << query_string;

     SqlQueryBuilder qbld;

     string sql_query = qbld.FormQuery(query_string);

     SLOG("NuVld", pINFO) << "Formed SQL-Query: " << sql_query;

     TSQLResult * res = fSqlServer->Query( sql_query.c_str() );

     const int nrows = res->GetRowCount();
     SLOG("NuVld", pINFO) << "The Query returned: " << nrows << " rows";

     if(nrows == 0) return eDbu_EMPTY_TABLE;
     
     TSQLRow *      row      = 0;
     DBSFTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new DBSFTableRow(row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;

     this  -> AddDBMeasurementIdList (table);
     table -> SetQueryString       (query_string);
     
     return eDbu_OK;     
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
void DBI::AddDBMeasurementIdList(DBTableBase * table)
{
// Builds a list of all the different measurements (i.e. different (experiment
// name, measurement tag) pairs) that can be found at this xsec_table.
//
  SLOG("NuVld", pINFO) << "Adding DBMeasurementIdList";

  DBMeasurementIdList * id_list = new DBMeasurementIdList();

  for(int irow = 0; irow < table->NRows(); irow++) {

     const DBTableRow * dbt_row = table->Row(irow);
  
     string experiment      = dbt_row->Field("name");
     string measurement_tag = dbt_row->Field("measurement_tag");

     DBMeasurementId * id = new DBMeasurementId(experiment, measurement_tag);

     ostringstream query;

     query << "SELECT authors, journal, year FROM REFERENCE "
           << "WHERE name = \"" << experiment.c_str() << "\""
           << " AND measurement_tag = \"" << measurement_tag.c_str() << "\"";

     SLOG("NuVld", pDEBUG) << "Query = " << query.str().c_str();

     TSQLResult * res = fSqlServer->Query(query.str().c_str());
     TSQLRow *    row = 0;

     unsigned int nrows = res->GetRowCount();

     for(unsigned int irow = 0; irow < nrows; irow++) {

         row = res->Next();

         id->_author.push_back  (string( row->GetField(0) ));
         id->_journal.push_back (string( row->GetField(1) ));
         id->_year.push_back    (string( row->GetField(2) ));

         SLOG("NuVld", pDEBUG) << row->GetField(0) << " "
                               << row->GetField(1) << " " << row->GetField(2);
         delete row;
     }
     delete res;

     id_list->AddId(id);
  }

  SLOG("NuVld", pDEBUG) << *id_list;

  table->SetDBMeasurementIdList(id_list);
}
//____________________________________________________________________________
bool DBI::HaveConnection(void) const
{
  return true;
}
//______________________________________________________________________________
