//_____________________________________________________________________________
/*!

\class    genie::nuvld::MultiLineMsgBox

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_MULTI_LINE_MSG_BOX_H_
#define _NUVLD_MULTI_LINE_MSG_BOX_H_

#include <string>
#include <vector>
#include <TApplication.h>
#include <TList.h>
#include <TVirtualX.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <RQ_OBJECT.h>

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

   void close_window (void) { delete this;               }
   void ok           (void) { _main->SendCloseMessage(); }

private:

   void build_multi_line_label      (const vector<string> * text);
   void position_relative_to_parent (const TGWindow * main);

   TGTransientFrame *  _main;
   TList *             _labels;
   TGButton *          _ok_button;
   TGLayoutHints *     _layout;

   ClassDef(MultiLineMsgBox, 0)
};

} // nuvld namespace
} // genie namespace

#endif


