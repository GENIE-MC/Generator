//_____________________________________________________________________________
/*!

\class    genie::nuvld::SysLogSingleton

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

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

class SysLogSingleton : public TObject {

public:

   friend class NuVldMainFrame;

   static SysLogSingleton * Instance();

   TGStatusBar *     StatusBar   (void) { return fStatusBar;   }
   TGHProgressBar *  ProgressBar (void) { return fProgressBar; }
   TGTextEdit *      Log         (void) { return fLog;         }
   
private:

   SysLogSingleton();
   SysLogSingleton(const SysLogSingleton & syslog);
   ~SysLogSingleton();

   static SysLogSingleton * fSelf;
   
   TGStatusBar *     fStatusBar;
   TGHProgressBar *  fProgressBar;
   TGTextEdit *      fLog;

   ClassDef(SysLogSingleton, 0)
};

} // nuvld namespace
} // genie namespace

#endif

