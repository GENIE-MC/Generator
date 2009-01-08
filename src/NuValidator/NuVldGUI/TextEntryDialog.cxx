//_____________________________________________________________________________
/*!

\class    genie::nuvld::TextEntryDialog

\brief    A simple text entry pop-up dialog

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/TextEntryDialog.h"

using namespace genie::nuvld;

ClassImp(TextEntryDialog)

//______________________________________________________________________________
TextEntryDialog::TextEntryDialog(const TGWindow * p, const TGWindow * main,
               UInt_t w, UInt_t h, UInt_t options, string title, string & txt) :
_txt(&txt)
{
  const char * cname = "genie::nuvld::TextEntryDialog";
  
  _main = new TGTransientFrame(p, main, w, h, options);
  
  _main->Connect("CloseWindow()", cname, this, "CloseWindow()");

  UInt_t lhints1 = kLHintsTop    | kLHintsCenterX | kLHintsExpandX;
  UInt_t lhints2 = kLHintsBottom | kLHintsLeft;
  UInt_t lhints3 = kLHintsBottom | kLHintsRight;
  
  _layout_1 = new TGLayoutHints(lhints1, 2, 2, 2, 2);
  _layout_2 = new TGLayoutHints(lhints2, 2, 2, 2, 2);
  _layout_3 = new TGLayoutHints(lhints3, 2, 2, 2, 2);

  _text_entry = new TGTextEntry(_main, new TGTextBuffer(35));

  _main->AddFrame (_text_entry, _layout_1);

  _buttons = new TGCompositeFrame(_main, 3, 3, kHorizontalFrame);

  _ok_button = new TGTextButton(_buttons, "&Ok", 1);
  _ok_button->Connect("Clicked()", cname, this, "Ok()");

  _cancel_button = new TGTextButton(_buttons, "&Cancel", 2);
  _cancel_button->Connect("Clicked()", cname, this, "Cancel()");

  _buttons->AddFrame (_ok_button,     _layout_2);
  _buttons->AddFrame (_cancel_button, _layout_3);

  _main->AddFrame (_buttons);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // pos. relative to the parent's window

  _main->SetWindowName(title.c_str());

  _main->MapWindow();
  
  gClient->WaitFor(_main);
}
//______________________________________________________________________________
TextEntryDialog::~TextEntryDialog()
{
   delete _ok_button;
   delete _cancel_button;
   delete _buttons;
   delete _text_entry;
   delete _layout_1;
   delete _layout_2;
   delete _layout_3;
   delete _main;
}
//______________________________________________________________________________
void TextEntryDialog::PositionRelativeToParent(const TGWindow * main)
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
void TextEntryDialog::Ok(void)
{
  *_txt =  string(_text_entry->GetBuffer()->GetString());

  _main->SendCloseMessage();
}
//______________________________________________________________________________
