//_____________________________________________________________________________
/*!

\class    genie::nuvld::MsgBox

\brief    A GUI Message Box

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _MSG_BOX_H_
#define _MSG_BOX_H_

#include <RQ_OBJECT.h>

class TGFrame;
class TGLabel;
class TGButton;

namespace genie {
namespace nuvld {

class MsgBox {

RQ_OBJECT("MsgBox")

public:
   MsgBox(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                        UInt_t options = kVerticalFrame, const char * txt = "");
   virtual ~MsgBox();

   void CloseWindow (void) { delete this;               }
   void Ok          (void) { fMain->SendCloseMessage(); }

private:

   void position_relative_to_parent(const TGWindow * main);

   TGTransientFrame * fMain;
   TGLabel *          fLabel;
   TGButton *         fOkBtn;
   TGLayoutHints *    fLayout;

   ClassDef(MsgBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

