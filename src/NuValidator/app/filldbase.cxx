//____________________________________________________________________________
/*!

\program  filldbase

\brief    Program to parse an XML data file and upload it to a Data Base

\synopsis filldbase -f filename.xml -h host -d dbase -u username -p passwd

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created October 05, 2004
*/
//____________________________________________________________________________

#include <sstream>
#include <string>

#include "TSQLServer.h"

#include "DBUtils/DBI.h"
#include "DBUtils/DBStatus.h"
#include "Messenger/Messenger.h"
#include "XmlParser/NuVldXmlParser.h"
#include "XmlParser/ParserStatus.h"

using std::ostringstream;
using std::string;

using genie::Messenger;
using namespace genie::nuvld;

const char * get_argument(int argc, char ** argv, const char * marker);

//__________________________________________________________________________________________
int main(int argc, char ** argv)
{
 Messenger * msg = Messenger::Instance();

 msg->SetPriorityLevel("NuVld", pDEBUG);

 const char * filename = get_argument(argc, argv, "-f");
 const char * hostname = get_argument(argc, argv, "-h");
 const char * dbase    = get_argument(argc, argv, "-d");
 const char * user     = get_argument(argc, argv, "-u");
 const char * password = get_argument(argc, argv, "-p");
 
 SLOG("filldbase", pINFO) << "parsing document.........: " << filename;
 SLOG("filldbase", pINFO) << "hostname.................: " << hostname;
 SLOG("filldbase", pINFO) << "dbase....................: " << dbase;
 SLOG("filldbase", pINFO) << "user.....................: " << user;
 SLOG("filldbase", pINFO) << "password.................: " << password;

 NuVldXmlParser xml_parser;

 xml_parser.ParseXmlDocument( filename );

 if(xml_parser.GetXmlParsingStatus() == eXml_OK) {

    SLOG("filldbase", pINFO) << "*** XML document successfully parsed";

 } else {

    SLOG("filldbase", pERROR) 
          << "Problems parsing XML document: " 
              << ParserStatus::AsString( xml_parser.GetXmlParsingStatus() );

    exit(1);
 }

 ostringstream db_url;
 db_url << "mysql://" << hostname << "/" << dbase;

 TSQLServer * sql_server = TSQLServer::Connect(db_url.str().c_str(), user, password);

 if( sql_server->IsConnected() ) {

    SLOG("filldbase", pINFO)
            << "Connected to RDBMS: " << sql_server->GetDBMS() 
            << " at host: "           << sql_server->GetHost() 
            << " / port: "            << sql_server->GetPort();
 } else {

    SLOG("filldbase", pERROR) << "Connection to RDBMS server failed";
 }

 DBI dbi(sql_server);
 DBStatus_t db_status = dbi.UploadXML( xml_parser.GetDataSet() );

 if(db_status != eDbu_OK) {

    SLOG("filldbase", pERROR) << "Failed to upload data into the dbase";
    exit(3);
 }

 return 0;
}
//__________________________________________________________________________________________
const char * get_argument(int argc, char ** argv, const char * marker)
{
  for(int iarg = 0; iarg < argc-1; iarg++) {
     string argument(argv[iarg]);
     if( argument.compare( marker ) == 0 ) return argv[++iarg];
  }
  return "";
}
//__________________________________________________________________________________________
