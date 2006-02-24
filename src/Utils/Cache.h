//____________________________________________________________________________
/*!

\class    genie::Cache

\brief    GENIE Cache Memory

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 26, 2004

*/
//____________________________________________________________________________

#ifndef _CACHE_H_
#define _CACHE_H_

#include <map>
#include <string>
#include <ostream>

#include <TFile.h>

using std::map;
using std::string;
using std::ostream;

namespace genie {

class CacheBranchI;

class Cache
{
public:

  static Cache * Instance(void);

  CacheBranchI * FindCacheBranch (string key);
  void           AddCacheBranch  (string key, CacheBranchI * branch);
  string         CacheBranchKey  (string k0, string k1="", string k2="") const;

  //-- print cache buffers
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Cache & cache);

private:

  //-- auto-load/save
  void AutoLoad      (void);
  void AutoSave      (void);
  void OpenCacheFile (void);

  //-- singleton instance
  static Cache * fInstance;

  //-- map of cache buffers & cache file
  map<string, CacheBranchI * > * fCacheMap;
  TFile *                        fCacheFile;

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
