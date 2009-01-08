//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiBrowserSingleton

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

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

class GuiBrowserSingleton : public TObject {

public:

   friend class NuVldMainFrame;

   static GuiBrowserSingleton * Instance();

   TRootEmbeddedCanvas * PlotBrowser(void) { return _e_canvas;   }
   TGTextEdit *          TextBrowser(void) { return _text_edit;  }

private:

   GuiBrowserSingleton();
   GuiBrowserSingleton(const GuiBrowserSingleton & browser);
   ~GuiBrowserSingleton();

   static GuiBrowserSingleton * _self;

   TRootEmbeddedCanvas * _e_canvas;
   TGTextEdit *          _text_edit;

   ClassDef(GuiBrowserSingleton, 0)
};

} // nuvld namespace
} // genie namespace

#endif

