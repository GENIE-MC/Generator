//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - January 12, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TList.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/GUIUtils.h"

using std::ostringstream;

//____________________________________________________________________________

void genie::utils::gui::FillListBox(TGListBox * lb, const char * lbitems[] )
{
  int i = 0;
  while( lbitems[i] )
  {
      lb->AddEntry(lbitems[i], i);
      i++;
  }
}
//____________________________________________________________________________
void genie::utils::gui::FillListBox(
                               TGListBox * lb, const vector<string> * lbitems)
{
  int i = 0;
  vector<string>::const_iterator lbiter;

  for(lbiter = lbitems->begin(); lbiter != lbitems->end(); ++lbiter) {

    if( lbiter->size() > 0 ) {

      lb->AddEntry( lbiter->c_str(), i);
      i++;
    }
  }
}
//____________________________________________________________________________
void genie::utils::gui::SelectAllListBoxEntries(TGListBox * lb)
{
  SLOG("GuiUtils", pDEBUG) << "Selecting all listbox entries";
  SLOG("GuiUtils", pDEBUG) << "n-entries = " << lb->GetNumberOfEntries();

  for(int i = 0; i < lb->GetNumberOfEntries(); i++) lb->Select(i);

  lb->SelectionChanged();
}
//____________________________________________________________________________
void genie::utils::gui::ResetAllListBoxSelections(TGListBox * lb)
{
  SLOG("GuiUtils", pDEBUG) << "Reseting all listbox entries";
  SLOG("GuiUtils", pDEBUG) << "n-entries = " << lb->GetNumberOfEntries();

  TList *     selected_entries = new TList();
  TGLBEntry * selected_entry   = 0;

  lb->GetSelectedEntries(selected_entries);

  TIter lbentry(selected_entries);

  while( (selected_entry = (TGLBEntry *) lbentry.Next()) )
                              lb->Select( selected_entry->EntryId() , kFALSE);

  delete selected_entries;

  lb->SelectionChanged();
}
//____________________________________________________________________________
string genie::utils::gui::ListBoxSelectionAsString(
                                       TGListBox * lb, const char * lbitems[])
{
  TList selections;

  lb->GetSelectedEntries( &selections  );

  TGLBEntry * lbentry  = 0;

  TIter selection_iter(&selections);

  int i = 0;

  ostringstream str_select;

  while( (lbentry = (TGLBEntry *) selection_iter.Next()) ) {

       str_select << lbitems[ lbentry->EntryId() ];

       if(++i < selections.GetSize() ) str_select << ", ";
  }

  if(i==0) return "empty";

  return str_select.str();
}
//____________________________________________________________________________
int genie::utils::gui::ListBoxSelectionId(
                               const char * lbitems[], const char * selection)
{
  int i = 0;
  while( lbitems[i] )
  {
      if ( strcmp(lbitems[i], selection) == 0 ) return i;
      i++;
  }
  return 0;
}
//____________________________________________________________________________
void genie::utils::gui::FillComboBox(TGComboBox * cb, const char * cbitems[])
{
  int i = 0;
  while( cbitems[i] )
  {
      cb->AddEntry(cbitems[i], i);
      i++;
  }
}
//____________________________________________________________________________
void genie::utils::gui::FillComboBox(
                              TGComboBox * cb, const vector<string> * cbitems)
{
  int i = 0;
  vector<string>::const_iterator cbiter;

  for(cbiter = cbitems->begin(); cbiter != cbitems->end(); ++cbiter) {

    if( cbiter->size() > 0 ) {
      cb->AddEntry( cbiter->c_str(), i);
      i++;
    }
  }
}
//____________________________________________________________________________
string genie::utils::gui::ComboBoxSelectionAsString(
                                     TGComboBox * cb, const char * cbitems[])
{
  TGLBEntry * selected_entry = cb->GetSelectedEntry();

  if(selected_entry) return string(cbitems[selected_entry->EntryId()]);
  else return "";
}
//____________________________________________________________________________
int genie::utils::gui::ComboBoxSelectionId(
                              const char * cbitems[], const char * selection)
{
  int i = 0;
  while( cbitems[i] ) {
  if ( strcmp(cbitems[i], selection) == 0 ) return i;
      i++;
  }
  return 0;
}
//____________________________________________________________________________


