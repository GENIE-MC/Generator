//_____________________________________________________________________________
/*!

\class    genie::nuvld::MsgBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _MSG_BOX_H_
#define _MSG_BOX_H_

#include <TApplication.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

namespace genie {
namespace nuvld {

class MsgBox {

RQ_OBJECT("MsgBox")

public:
   MsgBox(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
                        UInt_t options = kVerticalFrame, const char * txt = "");
   virtual ~MsgBox();

   void close_window (void) { delete this;               }
   void ok           (void) { _main->SendCloseMessage(); }

private:

   void position_relative_to_parent(const TGWindow * main);

   TGTransientFrame * _main;
   TGLabel *          _label;
   TGButton *         _ok_button;
   TGLayoutHints *    _layout;

   ClassDef(MsgBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif

