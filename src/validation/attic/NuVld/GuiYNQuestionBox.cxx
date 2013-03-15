//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TGLabel.h>
#include <TGButton.h>

#include "ValidationTools/NuVld/GuiYNQuestionBox.h"

using namespace genie::nuvld;

ClassImp(GuiYNQuestionBox)

//______________________________________________________________________________
GuiYNQuestionBox::GuiYNQuestionBox(const TGWindow *p,
             const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, 
                                               const char * txt, bool * answer):
fAnswer(answer)
{
  if(!txt) txt = "  ";

  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                          "genie::nuvld::GuiYNQuestionBox", this, "CloseWindow()");

  // add question label
  
  fLabel   = new TGLabel(fMain, new TGString(txt));  
  fLabelLt = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  // add control buttons
  
  fButtonFrame   = new TGHorizontalFrame(fMain, 10, 10);
  fButtonFrameLt = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 2,2,2,2);
  
  fYesBtn = new TGTextButton  (fButtonFrame, "&Yes", 1);
  fNoBtn  = new TGTextButton  (fButtonFrame, "&No",  2);
  
  fYesBtn -> Connect("Clicked()", "genie::nuvld::GuiYNQuestionBox", this, "Yes()");
  fNoBtn  -> Connect("Clicked()", "genie::nuvld::GuiYNQuestionBox", this, "No()" );

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
GuiYNQuestionBox::~GuiYNQuestionBox()
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
void GuiYNQuestionBox::PositionRelativeToParent(const TGWindow * main)
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
void GuiYNQuestionBox::Yes(void)
{
  *fAnswer = true;

  fMain->SendCloseMessage();
}
//______________________________________________________________________________
void GuiYNQuestionBox::No(void)
{
  *fAnswer = false;

  fMain->SendCloseMessage();
}
//______________________________________________________________________________

