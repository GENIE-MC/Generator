//____________________________________________________________________________
/*!

\class    genie::Cache

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 26, 2004

*/
//____________________________________________________________________________

#ifndef _CACHE_H_
#define _CACHE_H_

#include <map>
#include <string>

#include <TNtuple.h>

#include "Algorithm/Algorithm.h"

using std::map;
using std::string;

namespace genie {

class Cache
{
public:

  static Cache * Instance(void);

  string    CacheBranchKey     (const Algorithm * alg, string subr) const;
  TNtuple * FindCacheBranchPtr (const Algorithm * alg, string subr);
  TNtuple * CreateCacheBranch  (const Algorithm * alg, string subr, string brdef);
  
private:

  //-- singleton instance

  static Cache * fInstance;

  //-- map of cache buffers

  map<string, TNtuple * > * fCacheMap;

  //-- singleton class: constructors are private

  Cache();
  Cache(const Cache & cache);
  virtual ~Cache();

  //-- proper de-allocation of the singleton object

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (Cache::fInstance !=0) {
            delete Cache::fInstance;
            Cache::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _CACHE_H_
