//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiXmlFileHandler

\brief    Responds to GUI events associated with parsing, opening (through a
          TGFileDialog), closing, and uploading to an RDBMS of the NuValidator
          XML data files.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <vector>
#include <string>

#include <TSystem.h>
#include <TSQLResult.h>
#include <TGFileDialog.h>
#include <TGStatusBar.h>

#include "NuVldGUI/GuiXmlFileHandler.h"
#include "NuVldGUI/SysLogSingleton.h"
#include "NuVldGUI/MsgBox.h"
#include "NuVldGUI/MenuCommandId.h"

using std::vector;
using std::string;

using namespace genie::nuvld;

//______________________________________________________________________________
GuiXmlFileHandler::GuiXmlFileHandler()
{
  _main  = 0;
  _dbc   = 0;
}
//______________________________________________________________________________
GuiXmlFileHandler::GuiXmlFileHandler(const TGWindow * main, DBConnection * dbc):
_main(main),
_dbc(dbc)
{

}
//______________________________________________________________________________
GuiXmlFileHandler::~GuiXmlFileHandler()
{

}
//______________________________________________________________________________
void GuiXmlFileHandler::OpenFile(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  static TString dir(".");

  const char * k_xml_file_extensions[] =
  {
   "All files",     "*",
   "XML files",     "*.xml",
    0,               0
  };
  
  TGFileInfo fi;
  fi.fFileTypes = k_xml_file_extensions;
  fi.fIniDir    = StrDup( dir.Data() );

  new TGFileDialog(gClient->GetRoot(), _main, kFDOpen, &fi);

  if( fi.fFilename ) {

     _current_file = string( fi.fFilename );

     ostringstream cmd;
     cmd << "Opening file: " << _current_file.c_str();

     syslog -> Log()       -> AddLine ( cmd.str().c_str()    );
     syslog -> StatusBar() -> SetText ( cmd.str().c_str(), 0 );
     syslog -> StatusBar() -> SetText ( "XML File Open",   1 );

     _file_menu->EnableEntry(M_FILE_CLOSE);

     // fancy fake effect...
     for(int i=0; i<100; i++) {
//             gSystem->Sleep(3); syslog->progress_bar()->SetPosition(i);
     }
//     syslog->progress_bar()->SetPosition(0);
  }
}
//______________________________________________________________________________
void GuiXmlFileHandler::CloseFile(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  if( _current_file.size() == 0 ) {

     new MsgBox(gClient->GetRoot(), _main,
                                380, 250, kVerticalFrame, "No Input XML File");
  } else {

     ostringstream cmd;
     cmd << "Closing file: " << _current_file.c_str();

     syslog -> Log()       -> AddLine ( cmd.str().c_str() );
     syslog -> StatusBar() -> SetText( cmd.str().c_str(),   0 );
     syslog -> StatusBar() -> SetText( "No Input XML File", 1 );

     _current_file = "";

     _file_menu->DisableEntry(M_FILE_CLOSE);
  }
}
//______________________________________________________________________________
void GuiXmlFileHandler::ParseFile(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog -> Log()       -> AddLine( "Parsing XML data" );
  syslog -> StatusBar() -> SetText( "Parsing XML data", 0 );

  if( _current_file.size() == 0 ) {

      new MsgBox(gClient->GetRoot(), _main,
                                380, 250, kVerticalFrame, "No Input XML File");

      string mesg = "No Input XML File - Use File->Open first";

      syslog -> Log()       -> AddLine ( mesg.c_str()    );
      syslog -> StatusBar() -> SetText ( mesg.c_str(), 0 );

  } else {
    
      ostringstream cmd;
      cmd  << gSystem->Getenv("GENIE")
                               << "/bin/parsexml -f " << _current_file.c_str();
      ostringstream cmd2;
      cmd2  << "Executing: " << cmd.str().c_str();
      
      syslog -> Log()       ->AddLine ( cmd2.str().c_str() );
      syslog -> StatusBar() ->SetText ( cmd2.str().c_str(), 0 );

      // Run the 'parsexml' executable provided by the NuValidator package
      gSystem->Exec( cmd.str().c_str() );
  }
}
//______________________________________________________________________________
void GuiXmlFileHandler::SendToDBase(void)
{
  SysLogSingleton * syslog = SysLogSingleton::Instance();

  syslog -> Log()       -> AddLine( "Uploading parsed XML data to database"    );
  syslog -> StatusBar() -> SetText( "Uploading parsed XML data to database", 0 );

  if( ! this->IsConnected() )
         new MsgBox(gClient->GetRoot(), _main, 380, 250, kVerticalFrame,
                                "Undefinded DB - Use 'Connect to dbase' first");
  else {

     if(_current_file.size() == 0) {
          new MsgBox(gClient->GetRoot(), _main, 380, 250, kVerticalFrame,
                                  "No XML input file - Use 'File->Open' first");
     } else {

         ostringstream cmd;

         cmd  << gSystem->Getenv("GENIE") << "/bin/filldbase "
              << " -f " <<  _current_file.c_str()
              << " -h " <<  _dbc->Host()
              << " -d " <<  _dbc->DataBase()
              << " -u " <<  _dbc->User()
              << " -p " <<  _dbc->Password();

         ostringstream running_cmd;
         running_cmd << "Running: " << cmd.str().c_str();

         syslog -> Log()       -> AddLine( running_cmd.str().c_str() );
         syslog -> StatusBar() -> SetText( running_cmd.str().c_str(), 0 );

         // Run the filldbase stand-alone provided by the NuValidator package
         gSystem->Exec( cmd.str().c_str() );
     }
  }
}
//______________________________________________________________________________
bool GuiXmlFileHandler::IsConnected(void)
{
  if( !_dbc->SqlServer() ) return false;
  else 
    return  _dbc->SqlServer()->IsConnected();
}
//______________________________________________________________________________

