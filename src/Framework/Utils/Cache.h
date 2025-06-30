//____________________________________________________________________________
/*!

\class    genie::Cache

\brief    GENIE Cache Memory

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  November 26, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
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

class Cache;
class CacheBranchI;

ostream & operator << (ostream & stream, const Cache & cache);

class Cache
{
public:

  static Cache * Instance(void);

  //! cache file
  void OpenCacheFile (string filename);

  //! finding/adding cache branches
  CacheBranchI * FindCacheBranch (string key);
  void           AddCacheBranch  (string key, CacheBranchI * branch);
  string         CacheBranchKey  (string k0, string k1="", string k2="") const;

  //! removing cache branches
  void RmCacheBranch         (string key);
  void RmAllCacheBranches    (void);
  void RmMatchedCacheBranches(string key_substring);

  //! print cache buffers
  void   Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Cache & cache);

private:

  //! load/save
  void Load (void);
  void Save (void);

  //! singleton instance
  static Cache * fInstance;

  //! map of cache buffers & cache file
  map<string, CacheBranchI * > * fCacheMap;
  TFile *                        fCacheFile;

  //! singleton class: constructors are private
  Cache();
  Cache(const Cache & cache);
  virtual ~Cache();

  //! proper de-allocation of the singleton object
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
