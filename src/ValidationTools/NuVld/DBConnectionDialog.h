//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBConnectionDialog

\brief    A GUI dialog for connecting to an RDBMS

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _DBASE_CONNECTION_DIALOG_H_
#define _DBASE_CONNECTION_DIALOG_H_

#include <TGFrame.h>
#include <RQ_OBJECT.h>

class TGLabel;
class TGButton;
class TGTextEdit;
class TGTextEntry;
class TGGuiMsgBox;
class TGCheckButton;
class TGStatusBar;
class TGProgressBar;

namespace genie {
namespace nuvld {

class DBConnection;

class DBConnectionDialog {

RQ_OBJECT("DBConnectionDialog")

public:
   DBConnectionDialog(const TGWindow *p, const TGWindow *main, UInt_t w, 
                          UInt_t h, UInt_t options, DBConnection * connection);
   virtual ~DBConnectionDialog();

   void CloseWindow  (void)  { delete this; }
   void Cancel       (void);
   void OK           (void);

private:

   TGGroupFrame  *     BuildInputFieldsFrame          (void);
   TGHorizontalFrame * BuildControlButtonsFrame       (void);
   void                PositionRelativeToParent       (const TGWindow * main);
   void                CheckForDefaultsFromConfigFile (void);
   
   DBConnection *      _db_connection;

   TGTransientFrame *  _main;
   TGButton *          _ok_button;
   TGButton *          _cancel_button;
   TGTextEntry *       _url;
   TGTextEntry *       _user;
   TGTextEntry *       _dbase;
   TGTextEntry *       _passwd;
   TGCheckButton *     _keep_values;
   TGLayoutHints *     _text_layout;
   TGLayoutHints *     _button_layout;
   TGGroupFrame *      _grpf;
   TGHorizontalFrame * _button_frame;

   ClassDef(DBConnectionDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif

