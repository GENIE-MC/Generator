//_____________________________________________________________________________
/*!

\class    genie::nuvld::YNQuestionBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#include "NuVldGUI/YNQuestionBox.h"

using namespace genie::nuvld;

ClassImp(YNQuestionBox)

//______________________________________________________________________________
YNQuestionBox::YNQuestionBox(const TGWindow *p,
             const TGWindow *main, UInt_t w, UInt_t h, UInt_t options, 
                                               const char * txt, bool * answer):
_answer(answer)
{
  if(!txt) txt = "  ";

  _main = new TGTransientFrame(p, main, w, h, options);
  _main->Connect("CloseWindow()",
                          "genie::nuvld::YNQuestionBox", this, "CloseWindow()");

  // add question label
  
  _label        = new TGLabel(_main, new TGString(txt));  
  _label_layout = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 2, 2, 2, 2);

  // add control buttons
  
  _button_frame        = new TGHorizontalFrame(_main, 10, 10);
  _button_frame_layout = new TGLayoutHints(
                                  kLHintsBottom | kLHintsCenterX,  2, 2, 2, 2);
  
  _yes_button = new TGTextButton  (_button_frame, "&Yes", 1);
  _no_button  = new TGTextButton  (_button_frame, "&No",  2);
  
  _yes_button -> Connect("Clicked()",
                                  "genie::nuvld::YNQuestionBox", this, "Yes()");
  _no_button  -> Connect("Clicked()",
                                  "genie::nuvld::YNQuestionBox", this, "No()" );

  _button_frame -> AddFrame (_yes_button);
  _button_frame -> AddFrame (_no_button );

  // add label & button frame to "question box"
  
  _main -> AddFrame (_label,        _label_layout       );
  _main -> AddFrame (_button_frame, _button_frame_layout);
  
  _main->MapSubwindows();
  _main->Resize();

  this->PositionRelativeToParent(main); // position relative to the parent's window

  _main->SetWindowName("QuestionBox");

  _main->MapWindow();
  gClient->WaitFor(_main);
}
//______________________________________________________________________________
YNQuestionBox::~YNQuestionBox()
{
  delete _label_layout;
  delete _label;
  delete _yes_button;
  delete _no_button;
  delete _button_frame_layout;
  delete _button_frame;
  delete _main;
}
//______________________________________________________________________________
void YNQuestionBox::PositionRelativeToParent(const TGWindow * main)
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
void YNQuestionBox::Yes(void)
{
  *_answer = true;

  _main->SendCloseMessage();
}
//______________________________________________________________________________
void YNQuestionBox::No(void)
{
  *_answer = false;

  _main->SendCloseMessage();
}
//______________________________________________________________________________

