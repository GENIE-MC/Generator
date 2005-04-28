//_____________________________________________________________________________
/*!

\class    genie::nuvld::MultiLineMsgBox

\brief    A multi-line GUI Message Box

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <TGFrame.h>
#include <TGIcon.h>
#include <TList.h>
#include <TGLabel.h>
#include <TGButton.h>

#include "NuVldGUI/MultiLineMsgBox.h"

using namespace genie::nuvld;

ClassImp(MultiLineMsgBox)

//______________________________________________________________________________
MultiLineMsgBox::MultiLineMsgBox(const TGWindow *p,
                const TGWindow *main, UInt_t w, UInt_t h,
                                     UInt_t options, const vector<string> * txt)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                        "genie::nuvld::MultiLineMsgBox", this, "CloseWindow()");

  fLayout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  fLbl = new TList();
  
  if(txt) this->BuildMultiLineLabel(txt);   // add text lines
  
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::MultiLineMsgBox", this, "Ok()");

  fMain->AddFrame (fOkBtn,  fLayout);

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("MsgBox");

  fMain->MapWindow();
  gClient->WaitFor(fMain);
}
//______________________________________________________________________________
MultiLineMsgBox::~MultiLineMsgBox()
{
  fLbl->Clear();
  delete fLbl;
  delete fOkBtn;
  delete fLayout;
  delete fMain;
}
//______________________________________________________________________________
void MultiLineMsgBox::BuildMultiLineLabel(const vector<string> * txt)
{
  unsigned int i = 0;
  
  vector<string>::const_iterator txt_iter = txt->begin();

  for( ; txt_iter != txt->end(); ++txt_iter) {

    TGLabel * label = new TGLabel(fMain, new TGString( txt_iter->c_str() ));

    fLbl->Add(label);
    fMain->AddFrame (label, fLayout);
    i++;
  }
  
  fLbl->SetOwner(kTRUE);
}
//______________________________________________________________________________
void MultiLineMsgBox::PositionRelativeToParent(const TGWindow * main)
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

