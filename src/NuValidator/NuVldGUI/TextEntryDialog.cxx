//_____________________________________________________________________________
/*!

\class    genie::nuvld::TextEntryDialog

\brief    A simple text entry pop-up dialog

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/TextEntryDialog.h"

using namespace genie::nuvld;

ClassImp(TextEntryDialog)

//______________________________________________________________________________
TextEntryDialog::TextEntryDialog(
             const TGWindow * p,const TGWindow * main,
                         UInt_t w, UInt_t h, UInt_t options, const char * txt) :
_txt(txt)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()", "genie::nuvld::TextEntryDialog", this, "CloseWindow()");

  _layout_1 = new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 2, 2, 2);
  _layout_2 = new TGLayoutHints(kLHintsBottom | kLHintsLeft,   2, 2, 2, 2);
  _layout_3 = new TGLayoutHints(kLHintsBottom | kLHintsRight,  2, 2, 2, 2);

  _text_entry = new TGTextEntry(_main, new TGTextBuffer(350));

  _main->AddFrame (_text_entry, _layout_1);

  _buttons = new TGCompositeFrame(_main, 3, 3, kHorizontalFrame);

  _ok_button = new TGTextButton(_buttons, "&Ok", 1);
  _ok_button->Connect("Clicked()", "genie::nuvld::TextEntryDialog", this, "Ok()");

  _cancel_button = new TGTextButton(_buttons, "&Cancel", 2);
  _cancel_button->Connect("Clicked()", "genie::nuvld::TextEntryDialog", this, "Cancel()");

  _buttons->AddFrame (_ok_button,     _layout_2);
  _buttons->AddFrame (_cancel_button, _layout_3);

  _main->AddFrame (_buttons);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("Custom SQL Query");

  _main->MapWindow();
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
  _txt =  _text_entry -> GetBuffer() -> GetString();

  _main->SendCloseMessage();
}
//______________________________________________________________________________
