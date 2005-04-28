//_____________________________________________________________________________
/*!

\class    genie::nuvld::YNQuestionBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include <TGLabel.h>
#include <TGButton.h>

#include "NuVldGUI/YNQuestionBox.h"

using namespace genie::nuvld;

ClassImp(YNQuestionBox)

//______________________________________________________________________________
YNQuestionBox::YNQuestionBox(const TGWindow *p,
             const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, 
                                               const char * txt, bool * answer):
fAnswer(answer)
{
  if(!txt) txt = "  ";

  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                          "genie::nuvld::YNQuestionBox", this, "CloseWindow()");

  // add question label
  
  fLabel   = new TGLabel(fMain, new TGString(txt));  
  fLabelLt = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  // add control buttons
  
  fButtonFrame   = new TGHorizontalFrame(fMain, 10, 10);
  fButtonFrameLt = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 2,2,2,2);
  
  fYesBtn = new TGTextButton  (fButtonFrame, "&Yes", 1);
  fNoBtn  = new TGTextButton  (fButtonFrame, "&No",  2);
  
  fYesBtn -> Connect("Clicked()", "genie::nuvld::YNQuestionBox", this, "Yes()");
  fNoBtn  -> Connect("Clicked()", "genie::nuvld::YNQuestionBox", this, "No()" );

  fButtonFrame -> AddFrame (fYesBtn);
  fButtonFrame -> AddFrame (fNoBtn );

  // add label & button frame to "question box"
  
  fMain -> AddFrame (fLabel,        fLabelLt       );
  fMain -> AddFrame (fButtonFrame, fButtonFrameLt);
  
  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("QuestionBox");

  fMain->MapWindow();
  gClient->WaitFor(fMain);
}
//______________________________________________________________________________
YNQuestionBox::~YNQuestionBox()
{
  delete fLabelLt;
  delete fLabel;
  delete fYesBtn;
  delete fNoBtn;
  delete fButtonFrameLt;
  delete fButtonFrame;
  delete fMain;
}
//______________________________________________________________________________
void YNQuestionBox::PositionRelativeToParent(const TGWindow * main)
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
void YNQuestionBox::Yes(void)
{
  *fAnswer = true;

  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void YNQuestionBox::No(void)
{
  *fAnswer = false;

  fMain->SendCloseMessage();
}
//______________________________________________________________________________

