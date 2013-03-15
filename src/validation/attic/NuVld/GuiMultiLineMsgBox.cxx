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
#include <TGIcon.h>
#include <TList.h>
#include <TGLabel.h>
#include <TGButton.h>

#include "ValidationTools/NuVld/GuiMultiLineMsgBox.h"

using namespace genie::nuvld;

ClassImp(GuiMultiLineMsgBox)

//______________________________________________________________________________
GuiMultiLineMsgBox::GuiMultiLineMsgBox(const TGWindow *p,
                const TGWindow *main, UInt_t w, UInt_t h,
                                     UInt_t options, const vector<string> * txt)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()",
                 "genie::nuvld::GuiMultiLineMsgBox", this, "CloseWindow()");

  fLayout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  fLbl = new TList();
  
  if(txt) this->BuildMultiLineLabel(txt);   // add text lines
  
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::GuiMultiLineMsgBox", this, "Ok()");

  fMain->AddFrame (fOkBtn,  fLayout);

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("GuiMsgBox");

  fMain->MapWindow();
  gClient->WaitFor(fMain);
}
//______________________________________________________________________________
GuiMultiLineMsgBox::~GuiMultiLineMsgBox()
{
  fLbl->Clear();
  delete fLbl;
  delete fOkBtn;
  delete fLayout;
  delete fMain;
}
//______________________________________________________________________________
void GuiMultiLineMsgBox::BuildMultiLineLabel(const vector<string> * txt)
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
void GuiMultiLineMsgBox::PositionRelativeToParent(const TGWindow * main)
{
// position relative to the parent's window
  
  int ax, ay;
  Window_t wdum;

  gVirtualX->TranslateCoordinates(
             main->GetId(), fMain->GetParent()->GetId(),
             (Int_t)(((TGFrame *) main)->GetWidth()  - fMain->GetWidth())  >> 1,
             (Int_t)(((TGFrame *) main)->GetHeight() - fMain->GetHeight()) >> 1,
             ax, ay, wdum);

  fMain->Move(ax, ay);
}
//______________________________________________________________________________

