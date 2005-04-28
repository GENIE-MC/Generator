//_____________________________________________________________________________
/*!

\class    genie::nuvld::HelpBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <TGText.h>
#include <TGTextEdit.h>
#include <TGButton.h>
#include <TGWindow.h>

#include "NuVldGUI/HelpBox.h"

using namespace genie::nuvld;

ClassImp(HelpBox)

//______________________________________________________________________________
HelpBox::HelpBox(const TGWindow *p, const TGWindow *main,
                      UInt_t w, UInt_t h, UInt_t options, const char * filename)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()", "genie::nuvld::HelpBox", this, "CloseWindow()");

  fTextLt = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  fBtnLt = new TGLayoutHints( kLHintsTop | kLHintsCenterX, 2, 2, 2, 2 );
  fText  = new TGTextEdit(fMain,  400, 400, kSunkenFrame | kDoubleBorder);
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::HelpBox", this, "OK()");

  fMain->AddFrame (fText,       fTextLt  );
  fMain->AddFrame (fOkBtn,  fBtnLt);

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("HelpBox");

  this->LoadTextFromFle(filename);
  
  fMain->MapWindow();
  
  //gClient->WaitFor(fMain);
}
//______________________________________________________________________________
HelpBox::~HelpBox()
{
  delete fText;
  delete fOkBtn;
  delete fTextLt;
  delete fBtnLt;
  delete fMain;
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

    fText->AddText(&tgtxt);

    fclose(fp);
  }
}
//______________________________________________________________________________
void HelpBox::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window
  
  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(main->GetId(), fMain->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth() - fMain->GetWidth()) >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - fMain->GetHeight()) >> 1,
             ax, ay, wdum);
  fMain->Move(ax, ay);
}
//______________________________________________________________________________

