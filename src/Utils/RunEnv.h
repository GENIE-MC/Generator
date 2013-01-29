//____________________________________________________________________________
/*!

\class    genie::RunEnv

\brief    Run-time environment

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  January 29, 2013

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RUN_ENV_H_
#define _RUN_ENV_H_

#include <iostream>
#include <string>

using std::ostream;

namespace genie {

class RunEnv
{
public:
  static RunEnv * Instance(void);

  bool CacheEnabled(void) const { return fEnableCache; }

  void EnableCache(bool flag) { fEnableCache = flag; }

  // Print 
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const RunEnv & env);

private:

  //
  bool fEnableCache;

  // Self
  static RunEnv * fInstance;

  // Singleton class: constructors are private
  RunEnv();
  RunEnv(const RunEnv & env);
  virtual ~RunEnv();

  // Clean-up
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (RunEnv::fInstance !=0) {
            delete RunEnv::fInstance;
            RunEnv::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif // _RUN_ENV_H_
