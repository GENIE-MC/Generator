//____________________________________________________________________________
/*!

\class    genie::Cache

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 26, 2004

*/
//____________________________________________________________________________

#include <sstream>
#include <iostream>

#include "Messenger/Messenger.h"
#include "Utils/Cache.h"

using std::ostringstream;
using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
Cache * Cache::fInstance = 0;
//____________________________________________________________________________
Cache::Cache()
{
  fInstance =  0;
  fCacheMap = new map<string, TNtuple * >;
}
//____________________________________________________________________________
Cache::~Cache()
{
  cout << "Cache singleton dtor: Deleting all cache branches" << endl;
  fInstance = 0;
  map<string, TNtuple * >::iterator citer;
  for(citer = fCacheMap->begin(); citer != fCacheMap->end(); ++citer) {
    TNtuple * branch = citer->second;
    if(branch) {
      delete branch;
      branch = 0;
    }
  }
  fCacheMap->clear();
  delete fCacheMap;
}
//____________________________________________________________________________
Cache * Cache::Instance()
{
  if(fInstance == 0) {

    static Cache::Cleaner cleaner;

    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Cache;
  }
  return fInstance;
}
//____________________________________________________________________________
TNtuple * Cache::FindCacheBranchPtr(const Algorithm * alg, string subbranch)
{
// Each algorithm, in each of its configuration states, can have a number
// of cache branch to cache its data.
// The key for accessing the cache branch is assembled as:
// alg-name/alg-config-set/subbranch-number

  string key = this->CacheBranchKey(alg, subbranch);

  if (fCacheMap->count(key) == 1) {

     map<string, TNtuple *>::const_iterator map_iter;

     map_iter = fCacheMap->find(key);

     return map_iter->second;
  }

  return 0;
}
//____________________________________________________________________________
TNtuple * Cache::CreateCacheBranch(
                    const Algorithm * alg, string subbranch, string branchdef)
{
  string key = this->CacheBranchKey(alg, subbranch);

  TNtuple * nt = new TNtuple(key.c_str(), "cache branch", branchdef.c_str());

  nt->SetDirectory(0);
  nt->SetCircular(1600000);

  fCacheMap->insert( map<string, TNtuple *>::value_type(key,nt) );

  return nt;
}
//____________________________________________________________________________
string Cache::CacheBranchKey(const Algorithm * alg, string subbranch) const
{
  ostringstream key;

  key << alg->Name() << "/" << alg->ParamSet() << "/" << subbranch;

  return key.str();
}
//____________________________________________________________________________
