//_____________________________________________________________________________
/*!

\class    genie::nuvld::MenuCommandId

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _MENU_COMMAND_ID_H_
#define _MENU_COMMAND_ID_H_

namespace genie {
namespace nuvld {

enum EMenuCommandId {

   M_FILE_OPEN = 0,
   M_FILE_CLOSE,
   M_FILE_PARSE,
   M_LOAD_STACKED_DATA,
   M_SAVE_STACKED_DATA,  
   M_FILE_EXIT,
   M_DATA_CONNECT,
   M_DATA_CLOSE,
   M_DATA_VERIFY,
   M_DATA_DBINFO,
   M_DATA_BOOTSTRAP,
   M_DATA_UPLOAD,
   M_DATA_QUERY,
   M_DATA_QUERY_FILE,
   M_DATA_QUERY_DRAW_GUI,
   M_DATA_QUERY_PRINT_GUI,
   M_NEUGEN_CONFIG_PHYSICS,
   M_NEUGEN_CONFIG_PROCESS,
   M_NEUGEN_RUN,
   M_FIT_OPEN,
   M_FIT_RUN,
   M_FIT_RESET,
   M_SEND_SQL,
   M_SAVE_CANVAS,
   M_HELP_ABOUT,
   M_HELP_WWW,
   M_HELP_DURHAM,
   M_HELP_HOWTO_CONN_DBASE,
   M_HELP_HOWTO_FILL_DBASE
};

} // nuvld namespace
} // genie namespace

#endif

