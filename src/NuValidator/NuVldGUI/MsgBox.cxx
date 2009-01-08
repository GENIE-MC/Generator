//_____________________________________________________________________________
/*!

\class    genie::nuvld::MsgBox

\brief    A GUI Message Box

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>

#include "NuVldGUI/MsgBox.h"

using namespace genie::nuvld;

ClassImp(MsgBox)

//______________________________________________________________________________
MsgBox::MsgBox(const TGWindow *p,
     const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, const char * txt)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()", "genie::nuvld::MsgBox", this, "CloseWindow()");

  fLayout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  fLabel    = new TGLabel       (fMain, new TGString(txt));
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::MsgBox", this, "Ok()");

  fMain->AddFrame (fLabel,      fLayout);
  fMain->AddFrame (fOkBtn,  fLayout);

  fMain->MapSubwindows();
  fMain->Resize();

  position_relative_to_parent(main); // position relative to the parent's window

  fMain->SetWindowName("MsgBox");

  fMain->MapWindow();
  gClient->WaitFor(fMain);
}
//______________________________________________________________________________
MsgBox::~MsgBox()
{
  delete fLabel;
  delete fOkBtn;
  delete fLayout;
  delete fMain;
}
//______________________________________________________________________________
void MsgBox::position_relative_to_parent(const TGWindow * main)
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


