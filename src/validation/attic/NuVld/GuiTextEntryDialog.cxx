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

#include "ValidationTools/NuVld/GuiTextEntryDialog.h"

using namespace genie::nuvld;

ClassImp(GuiTextEntryDialog)

//______________________________________________________________________________
GuiTextEntryDialog::GuiTextEntryDialog(const TGWindow * p, const TGWindow * main,
               UInt_t w, UInt_t h, UInt_t options, string title, string & txt) :
_txt(&txt)
{
  const char * cname = "genie::nuvld::GuiTextEntryDialog";
  
  _main = new TGTransientFrame(p, main, w, h, options);
  
  _main->Connect("CloseWindow()", cname, this, "CloseWindow()");

  UInt_t lhints1 = kLHintsTop    | kLHintsCenterX | kLHintsExpandX;
  UInt_t lhints2 = kLHintsBottom | kLHintsLeft;
  UInt_t lhints3 = kLHintsBottom | kLHintsRight;
  
  _layout_1 = new TGLayoutHints(lhints1, 2, 2, 2, 2);
  _layout_2 = new TGLayoutHints(lhints2, 2, 2, 2, 2);
  _layout_3 = new TGLayoutHints(lhints3, 2, 2, 2, 2);

  _text_entry = new TGTextEntry(_main, new TGTextBuffer(35));

  _main->AddFrame (_text_entry, _layout_1);

  _buttons = new TGCompositeFrame(_main, 3, 3, kHorizontalFrame);

  _ok_button = new TGTextButton(_buttons, "&Ok", 1);
  _ok_button->Connect("Clicked()", cname, this, "Ok()");

  _cancel_button = new TGTextButton(_buttons, "&Cancel", 2);
  _cancel_button->Connect("Clicked()", cname, this, "Cancel()");

  _buttons->AddFrame (_ok_button,     _layout_2);
  _buttons->AddFrame (_cancel_button, _layout_3);

  _main->AddFrame (_buttons);

  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // pos. relative to the parent's window

  _main->SetWindowName(title.c_str());

  _main->MapWindow();
  
  gClient->WaitFor(_main);
}
//______________________________________________________________________________
GuiTextEntryDialog::~GuiTextEntryDialog()
{
   delete _ok_button;
   delete _cancel_button;
   delete _buttons;
   delete _text_entry;
   delete _layout_1;
   delete _layout_2;
   delete _layout_3;
   delete _main;
}
//______________________________________________________________________________
void GuiTextEntryDialog::PositionRelativeToParent(const TGWindow * main)
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
void GuiTextEntryDialog::Ok(void)
{
  *_txt =  string(_text_entry->GetBuffer()->GetString());

  _main->SendCloseMessage();
}
//______________________________________________________________________________
