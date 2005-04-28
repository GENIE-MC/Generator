//_____________________________________________________________________________
/*!

\class    genie::nuvld::MultiLineMsgBox

\brief    A multi-line GUI Message Box

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_MULTI_LINE_MSG_BOX_H_
#define _NUVLD_MULTI_LINE_MSG_BOX_H_

#include <string>
#include <vector>

#include <RQ_OBJECT.h>

class TGFrame;
class TGIcon;
class TList;
class TGLabel;
class TGButton;

using std::string;
using std::vector;

namespace genie {
namespace nuvld {

class MultiLineMsgBox {

RQ_OBJECT("MultiLineMsgBox")

public:
   MultiLineMsgBox(const TGWindow *p, const TGWindow *main, UInt_t w,
        UInt_t h, UInt_t options = kVerticalFrame, const vector<string> * text = 0);
   virtual ~MultiLineMsgBox();

   void CloseWindow (void)  { delete this;               }
   void Ok           (void) { fMain->SendCloseMessage(); }

private:

   void BuildMultiLineLabel      (const vector<string> * text);
   void PositionRelativeToParent (const TGWindow * main);

   TGTransientFrame *  fMain;
   TList *             fLbl;
   TGButton *          fOkBtn;
   TGLayoutHints *     fLayout;

   ClassDef(MultiLineMsgBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif


