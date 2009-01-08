//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiYNQuestionBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _Y_N_QUESTION_BOX_H_
#define _Y_N_QUESTION_BOX_H_

#include <TGFrame.h>
#include <RQ_OBJECT.h>

class TGLabel;
class TGButton;

namespace genie {
namespace nuvld {

class GuiYNQuestionBox {

RQ_OBJECT("GuiYNQuestionBox")

public:
   GuiYNQuestionBox(const TGWindow *p, const TGWindow *main, 
        UInt_t w, UInt_t h, UInt_t options, const char * txt, bool * answer);
   virtual ~GuiYNQuestionBox();

   void CloseWindow (void) { delete this; }
   void Yes         (void);
   void No          (void);

private:

   void PositionRelativeToParent(const TGWindow * main);

   bool *              fAnswer;
   TGTransientFrame *  fMain;
   TGHorizontalFrame * fButtonFrame;    
   TGLabel *           fLabel;
   TGButton *          fYesBtn;
   TGButton *          fNoBtn;
   TGLayoutHints *     fLabelLt;
   TGLayoutHints *     fButtonFrameLt;
   
   ClassDef(GuiYNQuestionBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

