//_____________________________________________________________________________
/*!

\class    genie::nuvld::MsgBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "MsgBox.h"

using namespace genie::nuvld;

ClassImp(MsgBox)

//namespace genie {
//namespace nuvld {

//______________________________________________________________________________
MsgBox::MsgBox(const TGWindow *p,
     const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, const char * txt)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()", "genie::nuvld::MsgBox", this, "close_window()");

  _layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  _label    = new TGLabel       (_main, new TGString(txt));
  _ok_button = new TGTextButton  (_main, "&Ok", 1);

  _ok_button->Connect("Clicked()", "genie::nuvld::MsgBox", this, "ok()");

  _main->AddFrame (_label,      _layout);
  _main->AddFrame (_ok_button,  _layout);

  _main->MapSubwindows();
  _main->Resize();

  position_relative_to_parent(main); // position relative to the parent's window

  _main->SetWindowName("MsgBox");

  _main->MapWindow();
  gClient->WaitFor(_main);
}
//______________________________________________________________________________
MsgBox::~MsgBox()
{
  delete _label;
  delete _ok_button;
  delete _layout;
  delete _main;
}
//______________________________________________________________________________
void MsgBox::position_relative_to_parent(const TGWindow * main)
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

//} // nuvld namespace
//} // genie namespace

