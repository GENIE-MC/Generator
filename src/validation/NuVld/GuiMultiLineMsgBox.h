//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiMultiLineMsgBox 

\brief    A multi-line GUI Message Box

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

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

class GuiMultiLineMsgBox  {

RQ_OBJECT("GuiMultiLineMsgBox")

public:
   GuiMultiLineMsgBox(const TGWindow *p, const TGWindow *main, UInt_t w,
        UInt_t h, UInt_t options = kVerticalFrame, const vector<string> * text = 0);
   virtual ~GuiMultiLineMsgBox();

   void CloseWindow (void)  { delete this;               }
   void Ok           (void) { fMain->SendCloseMessage(); }

private:

   void BuildMultiLineLabel      (const vector<string> * text);
   void PositionRelativeToParent (const TGWindow * main);

   TGTransientFrame *  fMain;
   TList *             fLbl;
   TGButton *          fOkBtn;
   TGLayoutHints *     fLayout;

   ClassDef(GuiMultiLineMsgBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif


