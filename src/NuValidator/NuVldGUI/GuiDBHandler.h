//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiDBHandler

\brief    Responds to GUI events associated with the NuValidator's data-base

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_DBASE_HANDLER_H_
#define _GUI_DBASE_HANDLER_H_

#include <string>

#include <TGWindow.h>
#include <TSQLResult.h>

#include "NuVldGUI/DBConnection.h"

using std::string;
using namespace genie::nuvld;

namespace genie {
namespace nuvld {

class GuiDBHandler {

public:

   GuiDBHandler();
   GuiDBHandler(const TGWindow * main, DBConnection * connection);
   ~GuiDBHandler();

   void MakeConnection         (void);
   void CloseConnection        (void);
   void CheckConnection        (void);
   void PrintInfo              (void);
   void Bootstrap              (void);
   void QueryWithSqlFromDialog (void);
   void QueryWithSqlFromFile   (void);
   bool IsConnected            (void);
   
private:

   string ReadSqlQueryFromFile       (string filename);
   void   PrintSqlResultInTGTextEdit (TSQLResult * rs);

   const TGWindow *   _main;
   DBConnection *     _dbc;
};

} // nuvld namespace
} // genie namespace

#endif

