//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiSysLogSingleton

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _SYSLOG_SINGLETON_H_
#define _SYSLOG_SINGLETON_H_

#include <TObject.h>

#include <TGStatusBar.h>
#include <TGTextEdit.h>

class TGHProgressBar;

namespace genie {
namespace nuvld {

class GuiSysLogSingleton : public TObject {

public:

   friend class NuVldMainFrame;

   static GuiSysLogSingleton * Instance();

   TGStatusBar *     StatusBar   (void) { return fStatusBar;   }
   TGHProgressBar *  ProgressBar (void) { return fProgressBar; }
   TGTextEdit *      Log         (void) { return fLog;         }
   
private:

   GuiSysLogSingleton();
   GuiSysLogSingleton(const GuiSysLogSingleton & syslog);
   ~GuiSysLogSingleton();

   static GuiSysLogSingleton * fSelf;
   
   TGStatusBar *     fStatusBar;
   TGHProgressBar *  fProgressBar;
   TGTextEdit *      fLog;

   ClassDef(GuiSysLogSingleton, 0)
};

} // nuvld namespace
} // genie namespace

#endif

