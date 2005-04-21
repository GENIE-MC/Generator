//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBI

\brief    Data Base Interface

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#include <TSQLResult.h>
#include <TSQLRow.h>

#include "DBUtils/DBI.h"
#include "DBUtils/DBXmlUploader.h"
#include "DBUtils/SqlQueryBuilder.h"
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
             DBTable<vXSecTableRow> * table, const DBQueryString & query_string)
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
     vXSecTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new vXSecTableRow(row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;
  
     this  -> AddMeasurementIdList (table);
     table -> SetQueryString       (query_string);

     return eDbu_OK;
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
DBStatus_t DBI::FillTable(
         DBTable<eDiffXSecTableRow> * table, const DBQueryString & query_string)
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
     
     TSQLRow *           row      = 0;
     eDiffXSecTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new eDiffXSecTableRow(row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;

     this  -> AddMeasurementIdList (table);
     table -> SetQueryString       (query_string);
     
     return eDbu_OK;     
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
DBStatus_t DBI::FillTable(
                DBTable<SFTableRow> * table, const DBQueryString & query_string)
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
     
     TSQLRow *    row      = 0;
     SFTableRow * xsec_row = 0;

     for (int i = 0; i < nrows; i++) {

       row      = res->Next();
       xsec_row = new SFTableRow(row);

       table->AddRow(xsec_row);

       delete row;
     }
     delete res;

     this  -> AddMeasurementIdList (table);
     table -> SetQueryString       (query_string);
     
     return eDbu_OK;     
  }

  SLOG("NuVld", pERROR) << "*** No connection to data-base";

  return eDbu_CONNECTION_LOST;
}
//______________________________________________________________________________
void DBI::AddMeasurementIdList(DBTableBase * table)
{
// Builds a list of all the different measurements (i.e. different (experiment
// name, measurement tag) pairs) that can be found at this xsec_table.
//
  SLOG("NuVld", pINFO) << "Adding MeasurementIdList";

  MeasurementIdList * id_list = new MeasurementIdList();

  for(int irow = 0; irow < table->NRows(); irow++) {

     const DBTableRow * dbt_row = table->Row(irow);
  
     string experiment      = dbt_row->Field("name");
     string measurement_tag = dbt_row->Field("measurement_tag");

     MeasurementId * id = new MeasurementId(experiment, measurement_tag);

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

  table->SetMeasurementIdList(id_list);
}
//____________________________________________________________________________
bool DBI::HaveConnection(void) const
{
  return true;
}
//______________________________________________________________________________
