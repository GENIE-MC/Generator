//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTextEntryDialog

\brief    A simple text entry pop-up dialog

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _TEXT_ENTRY_DIALOG_H_
#define _TEXT_ENTRY_DIALOG_H_

#include <iostream>
#include <string>
#include <sstream>
#include <TROOT.h>
#include <TApplication.h>
#include <TVirtualX.h>
#include <TSystem.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGTextEdit.h>
#include <TGStatusBar.h>
#include <RQ_OBJECT.h>

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;

namespace genie {
namespace nuvld {

class GuiTextEntryDialog {

RQ_OBJECT("GuiTextEntryDialog")

public:
   GuiTextEntryDialog(const TGWindow *p, const TGWindow *main,
        UInt_t w, UInt_t h, UInt_t options, string title, string & txt);
//        UInt_t w, UInt_t h, UInt_t options, const char * title, const char * txt);
   virtual ~GuiTextEntryDialog();

   void CloseWindow (void) { delete this;               }
   void Cancel      (void) { _main->SendCloseMessage(); }
   void Ok          (void);

private:

   void PositionRelativeToParent(const TGWindow * main);

//   const char * _txt;
   string * _txt;

   TGTransientFrame * _main;
   TGCompositeFrame * _buttons;
   TGButton *         _ok_button;
   TGButton *         _cancel_button;
   TGTextEntry *      _text_entry;
   TGLayoutHints *    _layout_1;
   TGLayoutHints *    _layout_2;
   TGLayoutHints *    _layout_3;

   ClassDef(GuiTextEntryDialog, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _TEXT_ENTRY_DIALOG_H_
