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

#include <TGText.h>
#include <TGTextEdit.h>
#include <TGButton.h>
#include <TGWindow.h>

#include "ValidationTools/NuVld/GuiHelpBox.h"

using namespace genie::nuvld;

ClassImp(GuiHelpBox)

//______________________________________________________________________________
GuiHelpBox::GuiHelpBox(const TGWindow *p, const TGWindow *main,
                      UInt_t w, UInt_t h, UInt_t options, const char * filename)
{
  fMain = new TGTransientFrame(p, main, w, h, options);
  fMain->Connect("CloseWindow()", "genie::nuvld::GuiHelpBox", this, "CloseWindow()");

  fTextLt = new TGLayoutHints(
     kLHintsTop | kLHintsCenterX | kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1);

  fBtnLt = new TGLayoutHints( kLHintsTop | kLHintsCenterX, 2, 2, 2, 2 );
  fText  = new TGTextEdit(fMain,  400, 400, kSunkenFrame | kDoubleBorder);
  fOkBtn = new TGTextButton  (fMain, "&Ok", 1);

  fOkBtn->Connect("Clicked()", "genie::nuvld::GuiHelpBox", this, "OK()");

  fMain->AddFrame (fText,       fTextLt  );
  fMain->AddFrame (fOkBtn,  fBtnLt);

  fMain->MapSubwindows();
  fMain->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  fMain->SetWindowName("GuiHelpBox");

  this->LoadTextFromFle(filename);
  
  fMain->MapWindow();
  
  //gClient->WaitFor(fMain);
}
//______________________________________________________________________________
GuiHelpBox::~GuiHelpBox()
{
  delete fText;
  delete fOkBtn;
  delete fTextLt;
  delete fBtnLt;
  delete fMain;
}
//______________________________________________________________________________
void GuiHelpBox::LoadTextFromFle(const char * filename)
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
void GuiHelpBox::PositionRelativeToParent(const TGWindow * main)
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

