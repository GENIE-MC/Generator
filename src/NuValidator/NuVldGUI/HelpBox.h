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

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGTextEdit.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

namespace genie {
namespace nuvld {

class HelpBox {

RQ_OBJECT("HelpBox")

public:
   HelpBox(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                        UInt_t options = kVerticalFrame, const char * txt = "");
   virtual ~HelpBox();

   void CloseWindow (void) { delete this;               }
   void OK          (void) { _main->SendCloseMessage(); }

private:

   void LoadTextFromFle          (const char * filename);
   void PositionRelativeToParent (const TGWindow * main);

   TGTransientFrame * _main;
   TGTextEdit *       _text;
   TGButton *         _ok_button;
   TGLayoutHints *    _text_layout;
   TGLayoutHints *    _button_layout;

   ClassDef(HelpBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

