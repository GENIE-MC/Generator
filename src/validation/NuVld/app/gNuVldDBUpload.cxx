//____________________________________________________________________________
/*!

\program  gnuvld_dbupload

\brief    Program to parse a NuValidator XML data file and upload it to a
          MySQL data-base.

\synopsis gnuvld_dbupload -f filename.xml -h host -d dbase -u user -p passwd

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created October 05, 2004
*/
//____________________________________________________________________________

#include <sstream>
#include <string>

#include <TSQLServer.h>

#include "Conventions/XmlParserStatus.h"
#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBStatus.h"
#include "ValidationTools/NuVld/NuVldXmlParser.h"

using std::ostringstream;
using std::string;

using namespace genie;
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
 
 SLOG("NuVld", pINFO) << "parsing document.........: " << filename;
 SLOG("NuVld", pINFO) << "hostname.................: " << hostname;
 SLOG("NuVld", pINFO) << "dbase....................: " << dbase;
 SLOG("NuVld", pINFO) << "user.....................: " << user;
 SLOG("NuVld", pINFO) << "password.................: " << password;

 NuVldXmlParser xml_parser;

 xml_parser.ParseXmlDocument( filename );

 if(xml_parser.GetXmlParsingStatus() == kXmlOK) {
    LOG("NuVld", pINFO) << "*** XML document successfully parsed";
 } else { 
    LOG("NuVld", pFATAL) 
        << "Problem! XML parser status : " 
        << XmlParserStatus::AsString( xml_parser.GetXmlParsingStatus() );
    exit(1);
 }

 ostringstream db_url;
 db_url << "mysql://" << hostname << "/" << dbase;

 TSQLServer * sql_server = TSQLServer::Connect(db_url.str().c_str(), user, password);

 if( sql_server->IsConnected() ) {
    SLOG("NuVld", pINFO)
            << "Connected to RDBMS: " << sql_server->GetDBMS() 
            << " at host: "           << sql_server->GetHost() 
            << " / port: "            << sql_server->GetPort();
 } else {
    SLOG("NuVld", pFATAL) 
            << "Connection to RDBMS server failed";
    exit(2);
 }

 DBI dbi(sql_server);
 DBStatus_t db_status = dbi.UploadXML( xml_parser.GetDataSet() );

 if(db_status != eDbu_OK) {
    SLOG("NuVld", pERROR) 
            << "Failed to upload data into the dbase";
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
