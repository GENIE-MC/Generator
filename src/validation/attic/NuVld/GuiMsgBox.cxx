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

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>

#include "ValidationTools/NuVld/GuiMsgBox.h"

using namespace genie::nuvld;

ClassImp(GuiMsgBox)

//______________________________________________________________________________
GuiMsgBox::GuiMsgBox(const TGWindow *p,
     const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, const char * txt)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()", "genie::nuvld::GuiMsgBox", this, "CloseWindow()");

  fLayout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  fLabel    = new TGLabel       (fMain, new TGString(txt));
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::GuiMsgBox", this, "Ok()");

  fMain->AddFrame (fLabel,      fLayout);
  fMain->AddFrame (fOkBtn,  fLayout);

  fMain->MapSubwindows();
  fMain->Resize();

  position_relative_to_parent(main); // position relative to the parent's window

  fMain->SetWindowName("GuiMsgBox");

  fMain->MapWindow();
  gClient->WaitFor(fMain);
}
//______________________________________________________________________________
GuiMsgBox::~GuiMsgBox()
{
  delete fLabel;
  delete fOkBtn;
  delete fLayout;
  delete fMain;
}
//______________________________________________________________________________
void GuiMsgBox::position_relative_to_parent(const TGWindow * main)
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


