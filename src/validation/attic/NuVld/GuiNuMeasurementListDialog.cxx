//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Nov 01, 2004

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

#include <cassert>
#include <sstream>
#include <vector>

#include <TGFrame.h>
#include <TGListBox.h>
#include <TGButton.h>
#include <TGText.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/DBConnection.h"
#include "ValidationTools/NuVld/GuiNuMeasurementListDialog.h"

using std::ostringstream;
using std::vector;

using namespace genie::nuvld;

ClassImp(GuiNuMeasurementListDialog)

//______________________________________________________________________________
GuiNuMeasurementListDialog::GuiNuMeasurementListDialog(
             const TGWindow *p, const TGWindow *main, bool * attn,
                         UInt_t w, UInt_t h, UInt_t options, DBConnection * db):
GuiDataSelectionDialog()
{
  _db = db;

  _attn  = attn;
  *_attn = true; // lock main window's attention through-out this dialog's lifetime
    
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                 "genie::nuvld::GuiNuMeasurementListDialog", this, "CloseWindow()");

  _listbox_layout = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  _button_layout = new TGLayoutHints( kLHintsTop | kLHintsCenterX,  2, 2, 2, 2);

  _measurements_listbox = new TGListBox    (_main,           201);
  _close_button         = new TGTextButton (_main, "&Close", 101);

  this->LoadXmlMeasurementsFromDB();

  _measurements_listbox -> Resize (580, 250);
  _measurements_listbox -> SetMultipleSelections( true );

  _close_button->Connect("Clicked()",
                       "genie::nuvld::GuiNuMeasurementListDialog", this, "Close()");

  _main->AddFrame (_measurements_listbox,  _listbox_layout  );
  _main->AddFrame (_close_button,          _button_layout);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("Cross Section XmlMeasurement Selection Dialog");

  _main->MapWindow();

  //gClient->WaitFor(_main);
}
//______________________________________________________________________________
GuiNuMeasurementListDialog::~GuiNuMeasurementListDialog()
{
  *_attn = false; // release attention-lock

  delete _measurements_listbox;
  delete _close_button;
  delete _listbox_layout;
  delete _button_layout;
  delete _main;
}
//______________________________________________________________________________
void GuiNuMeasurementListDialog::LoadXmlMeasurementsFromDB(void)
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

      LOG("NuVld", pINFO) << item.str();

      delete row;
  }
  delete result;
}
//______________________________________________________________________________
string GuiNuMeasurementListDialog::BundleKeyListInString(void)
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

      vector<string> key_elem = 
               utils::str::Split(entry->GetText()->GetString(),  ";");
      assert( key_elem.size() == 7 );

      key_list << key_elem[0] << "," << key_elem[1]; 
      if(ikey++ < nselected-1) key_list << ";";
  }

  return key_list.str().c_str();  
}
//______________________________________________________________________________
string GuiNuMeasurementListDialog::BundleCutsInString(void)
{
  return "";
}
//______________________________________________________________________________
string GuiNuMeasurementListDialog::BundleDrawOptInString(void)
{
  return "";
}
//______________________________________________________________________________
void GuiNuMeasurementListDialog::ResetSelections(void)
{

}
//______________________________________________________________________________
void GuiNuMeasurementListDialog::PositionRelativeToParent(const TGWindow * main)
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


