//____________________________________________________________________________
/*!

\class    genie::RunOpt

\brief    Some common run-time GENIE options.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  January 29, 2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RUN_OPT_H_
#define _RUN_OPT_H_

#include <iostream>
#include <string>

#include "Framework/Utils/TuneId.h"

class TBits;

using std::ostream;

namespace genie {

class RunOpt;
ostream & operator << (ostream & stream, const RunOpt & opt);

class RunOpt
{
public:
  static RunOpt * Instance(void);

  // Read options from the command line. Call from all GENIE command-line apps.
  void ReadFromCommandLine(int argc, char ** argv);

  // Get options set.
  TuneId * Tune                 (void) const { return fTune;                   }
  string EventGeneratorList     (void) const { return fEventGeneratorList;     }
  string CacheFile              (void) const { return fCacheFile;              }
  string MesgThresholdFiles     (void) const { return fMesgThresholds;         }
  TBits* UnphysEventMask        (void) const { return fUnphysEventMask;        }
  int    EventRecordPrintLevel  (void) const { return fEventRecordPrintLevel;  }
  int    MCJobStatusRefreshRate (void) const { return fMCJobStatusRefreshRate; }
  bool   BareXSecPreCalc        (void) const { return fEnableBareXSecPreCalc;  }
  string XMLPath                (void) const { return fXMLPath;  }

  // If a user accesses the GENIE objects directly, then most of the options above
  // can be set directly to the relevant objects (Messenger, Cache, etc).
  //
  void SetTuneName(string tuneName="Default");
  void BuildTune(); ///< build tune and inform XSecSplineList
  void SetEventGeneratorList(string evgenlist) { fEventGeneratorList = evgenlist; }
  void EnableBareXSecPreCalc(bool flag)        { fEnableBareXSecPreCalc = flag; }

  // Print
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const RunOpt & opt);

private:

  void Init (void);

  // options
  TuneId * fTune;                    ///< GENIE comprehensive neutrino interaction model tune.
  string fEventGeneratorList;        ///< Name of event generator list to be loaded by the event generation drivers.
  string fCacheFile;                 ///< Name of cache file, is cache is to be re-used.
  string fMesgThresholds;            ///< List of files (delimited with : if more than one) with custom mesg stream thresholds.
  TBits* fUnphysEventMask;           ///< Unphysical event mask.
  int    fEventRecordPrintLevel;     ///< GHEP event r ecord print level.
  int    fMCJobStatusRefreshRate;    ///< MC job status file refresh rate.
  bool   fEnableBareXSecPreCalc;     ///< Cache calcs relevant to free-nucleon xsecs before any nuclear xsec computation?
                                     ///< The option switches on/off cacheing calculations which interfere with event reweighting.
  string fXMLPath;                   ///< An path to look for XML in. Higher priority than GXMLPATH

  // Self
  static RunOpt * fInstance;

  // Singleton class: constructors are private
  RunOpt();
  RunOpt(const RunOpt & opt);
  virtual ~RunOpt();

  // Clean-up
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (RunOpt::fInstance !=0) {
            delete RunOpt::fInstance;
            RunOpt::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif // _RUN_OPT_H_
