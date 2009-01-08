//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiXmlFileHandler

\brief    Responds to GUI events associated with parsing, opening (through a
          TGFileDialog), closing, and uploading to an RDBMS of the NuValidator
          XML data files.

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_XML_FILE_HANDLER_H_
#define _GUI_XML_FILE_HANDLER_H_

class TGWindow;
class TGMenu;
class TGPopupMenu;

namespace genie {
namespace nuvld {

class GuiXmlFileHandler {

public:

   GuiXmlFileHandler();
   GuiXmlFileHandler(const TGWindow * main, DBConnection * connection);
   ~GuiXmlFileHandler();

   void OpenFile    (void);
   void CloseFile   (void);
   void ParseFile   (void);
   void SendToDBase (void);

   void AttachFileMenu(TGPopupMenu * file_menu) { _file_menu = file_menu; }
   
private:

   bool IsConnected (void);

   string             _current_file;
   TGPopupMenu *      _file_menu;
   const TGWindow *   _main;
   DBConnection *     _dbc;
};

} // nuvld namespace
} // genie namespace

#endif

