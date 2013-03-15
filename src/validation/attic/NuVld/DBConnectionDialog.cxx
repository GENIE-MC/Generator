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
*/
//____________________________________________________________________________ 

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <TSystem.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEdit.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGStatusBar.h>
#include <TGProgressBar.h>

#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/DBConnectionDialog.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "Utils/StringUtils.h"

using std::string;
using std::ostringstream;
using std::ios;
using std::ifstream;

using namespace genie::utils::str;
using namespace genie::nuvld;

ClassImp(DBConnectionDialog)

//______________________________________________________________________________
DBConnectionDialog::DBConnectionDialog(
    const TGWindow * p, const TGWindow * main, 
            UInt_t w, UInt_t h, UInt_t options, DBConnection * connection) :
_db_connection(connection)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                     "genie::nuvld::DBConnectionDialog", this, "CloseWindow()");

  //-- add text inputs

  _grpf        = this->BuildInputFieldsFrame();
  _text_layout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

  //-- add buttons

  _button_frame  = this->BuildControlButtonsFrame();
  _button_layout = new TGLayoutHints(kLHintsBottom | kLHintsCenterX,  2, 2, 2, 2);

  //-- add frames to main frame

  _main->AddFrame(_grpf,         _text_layout);
  _main->AddFrame(_button_frame, _button_layout);

  _main->MapSubwindows();
  _main->Resize();

  this->CheckForDefaultsFromConfigFile(); // defaults from user config file - if any
  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("Database Dialog");

  _main->MapWindow();

  gClient->WaitFor(_main);
}
//______________________________________________________________________________
DBConnectionDialog::~DBConnectionDialog()
{
  delete _url;
  delete _user;
  delete _dbase;
  delete _passwd;
  delete _text_layout;
  delete _grpf;
  delete _ok_button;
  delete _cancel_button;
  delete _button_layout;
  delete _button_frame;
  delete _main;
}
//______________________________________________________________________________
TGGroupFrame * DBConnectionDialog::BuildInputFieldsFrame(void)
{
  TGGroupFrame * gframe = new TGGroupFrame(_main, "Options", kVerticalFrame);
  
  gframe->SetTitlePos(TGGroupFrame::kRight); // right aligned
  gframe->SetLayoutManager(new TGMatrixLayout(gframe, 0, 2, 10));

  gframe->AddFrame(new TGLabel(gframe, new TGHotString("URL:")));
  gframe->AddFrame (_url    = new TGTextEntry(gframe, new TGTextBuffer(25)) );
  gframe->AddFrame(new TGLabel(gframe, new TGHotString("Data-base:")));
  gframe->AddFrame (_dbase  = new TGTextEntry(gframe, new TGTextBuffer(25)) );
  gframe->AddFrame(new TGLabel(gframe, new TGHotString("User:")));
  gframe->AddFrame (_user   = new TGTextEntry(gframe, new TGTextBuffer(25)) );
  gframe->AddFrame(new TGLabel(gframe, new TGHotString("Password:")));
  gframe->AddFrame (_passwd = new TGTextEntry(gframe, new TGTextBuffer(25)) );

  _url    -> Resize( 150, _url    -> GetDefaultHeight() );
  _dbase  -> Resize( 150, _dbase  -> GetDefaultHeight() );
  _user   -> Resize( 150, _user   -> GetDefaultHeight() );
  _passwd -> Resize( 150, _passwd -> GetDefaultHeight() );

  return gframe;
}
//______________________________________________________________________________
TGHorizontalFrame * DBConnectionDialog::BuildControlButtonsFrame(void)
{
  TGHorizontalFrame * bframe = new TGHorizontalFrame(_main, 10, 10);

  _ok_button = new TGTextButton(bframe, "&Ok", 1);
  _ok_button->Connect("Clicked()",
                              "genie::nuvld::DBConnectionDialog", this, "OK()");

  _cancel_button = new TGTextButton(bframe, "&Cancel", 2);
  _cancel_button->Connect("Clicked()",
                          "genie::nuvld::DBConnectionDialog", this, "Cancel()");

  bframe->AddFrame(_ok_button);
  bframe->AddFrame(_cancel_button);

  return bframe;
}
//______________________________________________________________________________
void DBConnectionDialog::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window

  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(main->GetId(), _main->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - _main->GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - _main->GetHeight()) >> 1,
             ax, ay, wdum);
  _main->Move(ax, ay);
}
//______________________________________________________________________________
void DBConnectionDialog::CheckForDefaultsFromConfigFile(void)
{
// set default values from user config file - if any

  string host, dbase, user, passwd;

  ifstream input(Concat(gSystem->Getenv("HOME"),"/.nuvld/dbase.vld"), ios::in);

  if(! input.fail() ) {

    input >> host;
    input >> dbase;
    input >> user;
    input >> passwd;

    input.close();
  }

  _url    -> SetText( host.c_str()   );
  _dbase  -> SetText( dbase.c_str()  );
  _user   -> SetText( user.c_str()   );
  _passwd -> SetText( passwd.c_str() );
}
//______________________________________________________________________________
void DBConnectionDialog::OK(void)
{
  GuiSysLogSingleton * syslog = GuiSysLogSingleton::Instance();

  string host    = string( _url    -> GetBuffer() -> GetString() );
  string dbase   = string( _dbase  -> GetBuffer() -> GetString() );
  string user    = string( _user   -> GetBuffer() -> GetString() );
  string passwd  = string( _passwd -> GetBuffer() -> GetString() );

  // display database information @ session log

  syslog->StatusBar()->SetText( "URL/Database/User/Password values entered", 0 );

  syslog->Log()->AddLine(Concat("URL........... ", host.c_str()   ) );
  syslog->Log()->AddLine(Concat("Database...... ", dbase.c_str()  ) );
  syslog->Log()->AddLine(Concat("User.......... ", user.c_str()   ) );
  syslog->Log()->AddLine(Concat("Passwd........ ", passwd.c_str() ) );

  // set DBConnection fields

  _db_connection->_host       = host;
  _db_connection->_database   = dbase;
  _db_connection->_user       = user;
  _db_connection->_password   = passwd;

  _main->SendCloseMessage();
}
//______________________________________________________________________________
void DBConnectionDialog::Cancel(void)
{
  _db_connection->_host       = "";
  _db_connection->_database   = "";
  _db_connection->_user       = "";
  _db_connection->_password   = "";
  _db_connection->_sql_server = 0;

  _main->SendCloseMessage();
}
//______________________________________________________________________________



