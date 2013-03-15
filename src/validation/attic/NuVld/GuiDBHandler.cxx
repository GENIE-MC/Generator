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

#include <cstdio>
#include <iostream>
#include <vector>

#include <TSystem.h>
#include <TGWindow.h>
#include <TGProgressBar.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TGFileDialog.h>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/GuiDBHandler.h"
#include "ValidationTools/NuVld/DBConnectionDialog.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiBrowserSingleton.h"
#include "ValidationTools/NuVld/GuiMultiLineMsgBox.h"
#include "ValidationTools/NuVld/GuiMsgBox.h"
#include "ValidationTools/NuVld/GuiYNQuestionBox.h"
#include "ValidationTools/NuVld/GuiTextEntryDialog.h"

using std::cout;
using std::endl;
using std::vector;

using namespace genie::utils::str;
using namespace genie::nuvld;

ClassImp(GuiDBHandler)

//______________________________________________________________________________
GuiDBHandler::GuiDBHandler()
{
  fMain = 0;
  fDBC  = 0;
}
//______________________________________________________________________________
GuiDBHandler::GuiDBHandler(const TGWindow * main, DBConnection * connection) : 
fMain(main),
fDBC(connection)
{

}
//______________________________________________________________________________
GuiDBHandler::~GuiDBHandler()
{

}
//______________________________________________________________________________
void GuiDBHandler::MakeConnection(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog -> Log()       ->AddLine( "Connecting to dbase" );
  syslog -> StatusBar() ->SetText( "Connecting to dbase", 0 );

  new DBConnectionDialog(
           gClient->GetRoot(), fMain, 800, 400, kVerticalFrame, fDBC);
            
  string url    = fDBC->URL();
  string user   = fDBC->User();
  string passwd = fDBC->Password();

  bool values_were_set = (url.size() > 0 && user.size() > 0);

  if(values_were_set) {
    
     TSQLServer * db = TSQLServer::Connect(
                                    url.c_str(), user.c_str(), passwd.c_str());
     if ( db ) {

        syslog -> Log() -> AddLine( "Server info: " );
        syslog -> Log() -> AddLine( db->ServerInfo() );

        fDBC->SetSqlServer(db);

     } else {

        syslog -> StatusBar() -> SetText("Connection to dbase failed", 0 );
        syslog -> Log()       -> AddLine("Connection to dbase failed");

        fDBC->SetSqlServer(0);

        new GuiMsgBox(gClient->GetRoot(), fMain,
                        380, 250, kVerticalFrame, "Connection to dbase failed");
     }
  }   
}
//______________________________________________________________________________
void GuiDBHandler::CloseConnection(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  if( ! this->IsConnected() )
    new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                            "There is no active dbase connection to terminate");
  else {

    // close connection

    syslog -> Log()       -> AddLine( "Closing dbase connection" );
    syslog -> StatusBar() -> SetText( "Closing dbase connection", 0 );
  
    fDBC->SqlServer()->Close();

    fDBC->SetSqlServer(0);
  }
}
//______________________________________________________________________________
void GuiDBHandler::CheckConnection(void)
{
  bool connected = this->IsConnected();

  if(connected)
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, 
                 kVerticalFrame, Concat(" Connected to ", fDBC->URL().c_str()));
  else
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                         "There is no active dbase connection");
}
//______________________________________________________________________________
void GuiDBHandler::PrintInfo(void)
{
  if(!this->IsConnected()) {

    new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                                        "There is no active dbase connection");
  } else {

    vector<string> db_info;

    TSQLServer * db = fDBC->SqlServer();

    db_info.push_back( "                                \
                                                 ");
    db_info.push_back( "   Getting SQL Server information:   " );
    db_info.push_back( "   -------------------------------   " );
    db_info.push_back( "     ");
    db_info.push_back( Concat("DBMS:..........", db->GetDBMS())    );
    db_info.push_back( Concat("Host:..........", db->GetHost())    );
    db_info.push_back( Concat("Port:..........", db->GetPort())    );
    db_info.push_back( Concat("Server info:...", db->ServerInfo()) );
    db_info.push_back( "     ");
    db_info.push_back( "  Tables:  " );
    db_info.push_back( "  -------  " );

    TSQLResult * tables = db->GetTables(fDBC->DataBase().c_str() );

    int nrows = tables->GetRowCount();

    for(int i=0; i<nrows; i++) {

      TSQLRow * table_row = tables->Next();

      TSQLResult * n_table_rows = db->Query( Concat(
                              "SELECT COUNT(*) FROM ", table_row->GetField(0)) );

      db_info.push_back(Concat(table_row->GetField(0), " ..... [nrows = ",
                                        n_table_rows->Next()->GetField(0),"]") );

      delete n_table_rows;
      delete table_row;
    }

    delete tables;

    db_info.push_back( "     ");

    new GuiMultiLineMsgBox(
           gClient->GetRoot(), fMain, 380, 250,  kVerticalFrame, &db_info);
  }
}
//______________________________________________________________________________
void GuiDBHandler::Bootstrap(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  if(! this->IsConnected()) {
      new GuiMsgBox(gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
                            "  Undefinded DB - Use 'Connect to dbase' first  ");
  } else {

     bool are_you_sure = false;

     // ask first if he really means to overwrite the database
     new GuiYNQuestionBox(
         gClient->GetRoot(), fMain, 380, 250, kVerticalFrame,
          "  I hope you were aware this will overwrite the dbase - Continue? ",
                                                                   &are_you_sure);

     // if the QuestionBox returns positive answer, bootstrap the dbase
     if(are_you_sure) {

        syslog -> Log()       -> AddLine( "Bootstraping the SQL data-base"    );
        syslog -> StatusBar() -> SetText( "Bootstraping the SQL data-base", 0 );

         fDBC->SqlServer()->DropDataBase("NuScat");

         fDBC->SqlServer()->CreateDataBase("NuScat");

         const int k_n_files = 7;

         string k_sql_file[k_n_files] = {
                            "createTable_ExpInfo.sql",
                            "createTable_BeamFlux.sql",
                            "createTable_Reference.sql",
                            "createTable_XmlMeasurementHeader.sql",
                            "createTable_CrossSection.sql",
                            "createTable_eDiffCrossSection.sql",
                            "createTable_StructureFunction.sql" };

         for(int ifile = 0; ifile < k_n_files; ifile++) {

            string filename = string(gSystem->Getenv("GENIE")) +
                                       "/data/sql_queries/" + k_sql_file[ifile];

            string sql = this->ReadSqlQueryFromFile(filename);
            LOG("NuVld", pINFO) << sql;
            
            fDBC->SqlServer()->Query( sql.c_str() );
         }
     }
  }
}
//______________________________________________________________________________
void GuiDBHandler::QueryWithSqlFromDialog(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();
  
  syslog -> Log()       -> AddLine( "Entering custom SQL query"    );
  syslog -> StatusBar() -> SetText( "Entering custom SQL query", 0 );

/*  
  new GuiTextEntryDialog(gClient->GetRoot(), fMain, 900, 500, sql);
  
  TSQLResult * result = fDBC->SqlServer()->Query( sql );

  syslog -> ProgressBar() -> SetPosition(0);
   
  this->PrintSqlResultInTGTextEdit(result);  

  syslog -> ProgressBar() -> SetPosition(0);
*/  
}
//______________________________________________________________________________
void GuiDBHandler::QueryWithSqlFromFile(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  syslog -> Log()       -> AddLine( "Loading custom SQL query from file" );
  syslog -> StatusBar() -> SetText( "Loading custom SQL query from file", 0 );

  static TString dir(".");

  const char * kSqlFileExt[] = {"All files", "*", "SQL files", "*.sql", 0, 0};
  
  TGFileInfo fi;
  fi.fFileTypes = kSqlFileExt;
  fi.fIniDir    = StrDup(dir.Data());

  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);

  if( fi.fFilename ) {

     string sqlFile = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Opening file: " << sqlFile.c_str();
     
     syslog -> Log()       -> AddLine( cmd.str().c_str() );
     syslog -> StatusBar() -> SetText( cmd.str().c_str(), 0 );
     syslog -> StatusBar() -> SetText( "SQL File Open", 1 );

     string sql = this->ReadSqlQueryFromFile(sqlFile);

     syslog -> ProgressBar() -> SetPosition(0);

     TSQLResult * result = fDBC->SqlServer()->Query( sql.c_str() );

     this->PrintSqlResultInTGTextEdit(result);     

     syslog -> ProgressBar() -> SetPosition(0);
  }
}
//______________________________________________________________________________
bool GuiDBHandler::IsConnected(void)
{
  if( !fDBC->SqlServer() ) return false;
  else 
    return fDBC->SqlServer()->IsConnected();
}
//______________________________________________________________________________
string GuiDBHandler::ReadSqlQueryFromFile(string filename)
{
  // read SQL query
  FILE *fp = fopen(filename.c_str(), "r");

  if(!fp) {
     //cerr << "File " << filename << " could not be read" << endl;
     return "";
  }

  char sql[4096] = {};
  fread(sql, 1, 4096, fp);
  fclose(fp);

  string ssql = string(sql);

  return utils::str::FilterString(";", ssql);
}
//______________________________________________________________________________
void GuiDBHandler::PrintSqlResultInTGTextEdit(TSQLResult * res)
{
  TSQLRow * row = 0;

  GuiSysLogSingleton * syslog   = GuiSysLogSingleton::Instance();
  GuiBrowserSingleton * browser = GuiBrowserSingleton::Instance();

  const int nr = res->GetRowCount();
  const int nf = res->GetFieldCount();

  string * field_name = new string[nf];

  for (int i = 0; i < nf; i++) field_name[i] = string( res->GetFieldName(i) );

  if(nr > 0) {
     double dprogress = 100. / nr;

     for (int i = 0; i < nr; i++) {

       syslog -> ProgressBar() -> SetPosition( (int) i*dprogress );

       row = res->Next();

       browser->TextBrowser()->AddLine(
                                  Concat("---------------------- row: ", i) );
       // print all fields
       for (int j = 0; j < nf; j++) {
           browser->TextBrowser()->AddLine(Concat(
                           field_name[j].c_str(), " : ", row->GetField(j) ) );
       }//fields
     }//rows
  }
  
  delete [] field_name;
}
//______________________________________________________________________________
