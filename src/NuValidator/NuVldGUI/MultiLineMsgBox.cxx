//_____________________________________________________________________________
/*!

\class    genie::nuvld::MultiLineMsgBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "MultiLineMsgBox.h"

using namespace genie::nuvld;

ClassImp(MultiLineMsgBox)

//namespace genie {
//namespace nuvld {

//______________________________________________________________________________
MultiLineMsgBox::MultiLineMsgBox(const TGWindow *p,
                const TGWindow *main, UInt_t w, UInt_t h,
                                     UInt_t options, const vector<string> * txt)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()", "genie::nuvld::MultiLineMsgBox", this, "close_window()");

  _layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  _labels = new TList();
  
  if(txt) build_multi_line_label(txt);   // add text lines
  
  _ok_button = new TGTextButton  (_main, "&Ok", 1);

  _ok_button->Connect("Clicked()", "genie::nuvld::MultiLineMsgBox", this, "ok()");

  _main->AddFrame (_ok_button,  _layout);

  _main->MapSubwindows();
  _main->Resize();

  position_relative_to_parent(main); // position relative to the parent's window

  _main->SetWindowName("MsgBox");

  _main->MapWindow();
  gClient->WaitFor(_main);
}
//______________________________________________________________________________
MultiLineMsgBox::~MultiLineMsgBox()
{
  _labels->Clear();
  delete _labels;
  delete _ok_button;
  delete _layout;
  delete _main;
}
//______________________________________________________________________________
void MultiLineMsgBox::build_multi_line_label(const vector<string> * txt)
{
  unsigned int i = 0;
  
  vector<string>::const_iterator txt_iter = txt->begin();

  for( ; txt_iter != txt->end(); ++txt_iter) {

    TGLabel * label = new TGLabel(_main, new TGString( txt_iter->c_str() ));

    _labels->Add(label);
    _main->AddFrame (label, _layout);
    i++;
  }
  
  _labels->SetOwner(kTRUE);
}
//______________________________________________________________________________
void MultiLineMsgBox::position_relative_to_parent(const TGWindow * main)
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
