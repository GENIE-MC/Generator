//_____________________________________________________________________________
/*!

\class    genie::nuvld::HelpBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <TGText.h>

#include "NuVldGUI/HelpBox.h"

using namespace genie::nuvld;

ClassImp(HelpBox)

//______________________________________________________________________________
HelpBox::HelpBox(const TGWindow *p, const TGWindow *main,
                      UInt_t w, UInt_t h, UInt_t options, const char * filename)
{
  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()", "genie::nuvld::HelpBox", this, "CloseWindow()");

  _text_layout = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  _button_layout = new TGLayoutHints( kLHintsTop | kLHintsCenterX, 2, 2, 2, 2 );

  _text      = new TGTextEdit(_main,  400, 400, kSunkenFrame | kDoubleBorder);
  _ok_button = new TGTextButton  (_main, "&Ok", 1);

  _ok_button->Connect("Clicked()", "genie::nuvld::HelpBox", this, "OK()");

  _main->AddFrame (_text,       _text_layout  );
  _main->AddFrame (_ok_button,  _button_layout);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("HelpBox");

  this->LoadTextFromFle(filename);
  
  _main->MapWindow();
  
  //gClient->WaitFor(_main);
}
//______________________________________________________________________________
HelpBox::~HelpBox()
{
  delete _text;
  delete _ok_button;
  delete _text_layout;
  delete _button_layout;
  delete _main;
}
//______________________________________________________________________________
void HelpBox::LoadTextFromFle(const char * filename)
{
  FILE * fp = fopen(filename, "r");

  char txt_buf[32000] = {};

  if(fp) {

    fread(txt_buf, 1, 32000, fp);

    TGText tgtxt;
    tgtxt.LoadBuffer(txt_buf);

    _text->AddText(&tgtxt);

    fclose(fp);
  }
}
//______________________________________________________________________________
void HelpBox::PositionRelativeToParent(const TGWindow * main)
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

