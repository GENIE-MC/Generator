//_____________________________________________________________________________
/*!

\class    genie::nuvld::BrowserSingleton

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created January 20, 2004
*/
//_____________________________________________________________________________

#ifndef _BROWSER_SINGLETON_H_
#define _BROWSER_SINGLETON_H_

#include <TObject.h>

class TRootEmbeddedCanvas;
class TGTextEdit;

namespace genie {
namespace nuvld {

class BrowserSingleton : public TObject {

public:

   friend class NuVldMainFrame;

   static BrowserSingleton * Instance();

   TRootEmbeddedCanvas * PlotBrowser(void) { return _e_canvas;   }
   TGTextEdit *          TextBrowser(void) { return _text_edit;  }

private:

   BrowserSingleton();
   BrowserSingleton(const BrowserSingleton & browser);
   ~BrowserSingleton();

   static BrowserSingleton * _self;

   TRootEmbeddedCanvas * _e_canvas;
   TGTextEdit *          _text_edit;

   ClassDef(BrowserSingleton, 0)
};

} // nuvld namespace
} // genie namespace

#endif

