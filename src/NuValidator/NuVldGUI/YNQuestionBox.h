//_____________________________________________________________________________
/*!

\class    genie::nuvld::YNQuestionBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _Y_N_QUESTION_BOX_H_
#define _Y_N_QUESTION_BOX_H_

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

namespace genie {
namespace nuvld {

class YNQuestionBox {

RQ_OBJECT("YNQuestionBox")

public:
   YNQuestionBox(const TGWindow *p, const TGWindow *main, 
        UInt_t w, UInt_t h, UInt_t options, const char * txt, bool * answer);
   virtual ~YNQuestionBox();

   void CloseWindow (void) { delete this; }
   void Yes         (void);
   void No          (void);

private:

   void PositionRelativeToParent(const TGWindow * main);

   bool *              _answer;

   TGTransientFrame *  _main;
   TGHorizontalFrame * _button_frame;    
   TGLabel *           _label;
   TGButton *          _yes_button;
   TGButton *          _no_button;
   TGLayoutHints *     _label_layout;
   TGLayoutHints *     _button_frame_layout;
   
   ClassDef(YNQuestionBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

