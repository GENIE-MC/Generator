//_____________________________________________________________________________
/*!

\class    genie::nuvld::HelpBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _HELP_BOX_H_
#define _HELP_BOX_H_

#include <TGFrame.h>
#include <RQ_OBJECT.h>

class TGTextEdit;
class TGButton;
class TGTransientFrame;

namespace genie {
namespace nuvld {

class HelpBox {

RQ_OBJECT("HelpBox")

public:
   HelpBox(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                        UInt_t options = kVerticalFrame, const char * txt = "");
   virtual ~HelpBox();

   void CloseWindow (void) { delete this;               }
   void OK          (void) { fMain->SendCloseMessage(); }

private:

   void LoadTextFromFle          (const char * filename);
   void PositionRelativeToParent (const TGWindow * main);

   TGTransientFrame * fMain;
   TGTextEdit *       fText;
   TGButton *         fOkBtn;
   TGLayoutHints *    fTextLt;
   TGLayoutHints *    fBtnLt;

   ClassDef(HelpBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

