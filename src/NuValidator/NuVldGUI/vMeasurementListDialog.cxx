//_____________________________________________________________________________
/*!

\class    genie::nuvld::vMeasurementListDialog

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  November 01, 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

#include <TGText.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "NuVldGUI/vMeasurementListDialog.h"
#include "XmlParser/ParserUtils.h"

using std::ostringstream;
using std::cout;
using std::endl;
using std::vector;
using std::string;

using namespace genie::nuvld;

ClassImp(vMeasurementListDialog)

//______________________________________________________________________________
vMeasurementListDialog::vMeasurementListDialog(
             const TGWindow *p, const TGWindow *main, bool & attn,
                         UInt_t w, UInt_t h, UInt_t options, DBConnection * db):
DataSelectionDialog(attn)
{
  _db = db;
  
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                 "genie::nuvld::vMeasurementListDialog", this, "CloseWindow()");

  _listbox_layout = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  _button_layout = new TGLayoutHints( kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);

  _measurements_listbox = new TGListBox    (_main,           201);
  _close_button         = new TGTextButton (_main, "&Close", 101);

  this->LoadMeasurementsFromDB();

  _measurements_listbox -> Resize (580, 250);
  _measurements_listbox -> SetMultipleSelections( true );

  _close_button->Connect("Clicked()",
                       "genie::nuvld::vMeasurementListDialog", this, "Close()");

  _main->AddFrame (_measurements_listbox,  _listbox_layout  );
  _main->AddFrame (_close_button,          _button_layout);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("Cross Section Measurement Selection Dialog");

  _main->MapWindow();

  //gClient->WaitFor(_main);
}
//______________________________________________________________________________
vMeasurementListDialog::~vMeasurementListDialog()
{
  delete _measurements_listbox;
  delete _close_button;
  delete _listbox_layout;
  delete _button_layout;
  delete _main;
}
//______________________________________________________________________________
void vMeasurementListDialog::LoadMeasurementsFromDB(void)
{
  TSQLServer * sql_server = _db->SqlServer();

  const char query[] =
        "SELECT MEASUREMENT_HEADER.observable, MEASUREMENT_HEADER.reaction, \
         MEASUREMENT_HEADER.name, MEASUREMENT_HEADER.measurement_tag, \
         REFERENCE.authors, REFERENCE.journal, REFERENCE.year \
         FROM MEASUREMENT_HEADER, REFERENCE \
         WHERE REFERENCE.name = MEASUREMENT_HEADER.name AND \
         REFERENCE.measurement_tag = MEASUREMENT_HEADER.measurement_tag \
         AND MEASUREMENT_HEADER.reaction LIKE \"%nu%\";";

  TSQLResult * result = sql_server->Query(query);

  const int nrows = result->GetRowCount();

  TSQLRow * row = 0;

  for (int i = 0; i < nrows; i++) {

      row = result->Next();

      ostringstream item;

      item << row->GetField(2) << ";"
           << row->GetField(3) << ";"
           << row->GetField(4) << ";"
           << row->GetField(5) << ";"
           << row->GetField(6) << ";"
           << row->GetField(0) << ";"
           << row->GetField(1) << ";";

      _measurements_listbox->AddEntry(item.str().c_str(), i);

      cout << item.str() << endl;

      delete row;
  }
  delete result;
}
//______________________________________________________________________________
string vMeasurementListDialog::BundleKeyListInString(void)
{
  int ikey = 0;
  
  ostringstream key_list;
  
  TList * selected = new TList();

  _measurements_listbox->GetSelectedEntries(selected);

  TGTextLBEntry * entry = 0;

  TIter iter(selected);

  // IndexOf() is broken in ROOT > 4.02 ??
  //int nselected = selected->IndexOf( selected->Last() ) + 1;
  int nselected = 0;
  while( (entry = (TGTextLBEntry *) iter.Next()) ) nselected++;
  iter.Reset();
  
  while( (entry = (TGTextLBEntry *) iter.Next()) ) {

      vector<string> key_elem = ParserUtils::split(
                                     entry->GetText()->GetString(),  ";");

      assert( key_elem.size() == 7 );

      key_list << key_elem[0] << "," << key_elem[1]; 

      if(ikey++ < nselected-1) key_list << ";";
  }

  return key_list.str().c_str();  
}
//______________________________________________________________________________
string vMeasurementListDialog::BundleCutsInString(void)
{
  return "";
}
//______________________________________________________________________________
string vMeasurementListDialog::BundleDrawOptInString(void)
{
  return "";
}
//______________________________________________________________________________
void vMeasurementListDialog::PositionRelativeToParent(const TGWindow * main)
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


