//____________________________________________________________________________
/*!

\namespace  genie::utils::gui

\brief      Simple utilities for GENIE Graphical User Interface widgets

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    January 12, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GUI_UTILS_H_
#define _GUI_UTILS_H_

#include <vector>
#include <string>

#include <TGListBox.h>
#include <TGComboBox.h>

using std::vector;
using std::string;

namespace genie {
namespace utils {
namespace gui  
{
  //-- ListBox methods

  void   FillListBox               (TGListBox * lb, const char * lbitems[]);
  void   FillListBox               (TGListBox * lb, const vector<string> * lbitems);
  void   SelectAllListBoxEntries   (TGListBox * lb);
  void   ResetAllListBoxSelections (TGListBox * lb);
  string ListBoxSelectionAsString  (TGListBox * lb, const char * lbitems[]);
  int    ListBoxSelectionId        (const char * lbitems[], const char * sel);

  //-- ComboBox methods

  void   FillComboBox              (TGComboBox * cb, const char * cbitems[]);
  void   FillComboBox              (TGComboBox * cb, const vector<string> * cbitems);
  string ComboBoxSelectionAsString (TGComboBox * cb, const char * cbitems[]);
  int    ComboBoxSelectionId       (const char * cbitems[], const char * sel);

}      // gui   namespace
}      // utils namespace
}      // genie namespace

#endif // _GUI_UTILS_H_
