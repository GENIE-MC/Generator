//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 26, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>
#include <iostream>

#include <TSystem.h>
#include <TDirectory.h>
#include <TList.h>
#include <TObjString.h>

#include "Messenger/Messenger.h"
#include "Utils/Cache.h"
#include "Utils/CacheBranchI.h"

using std::ostringstream;
using std::cout;
using std::endl;

namespace genie {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const Cache & cache)
{
  cache.Print(stream);
  return stream;
}
//____________________________________________________________________________
Cache * Cache::fInstance = 0;
//____________________________________________________________________________
Cache::Cache()
{
  fInstance  = 0;
  fCacheMap  = 0;
  fCacheFile = 0;
}
//____________________________________________________________________________
Cache::~Cache()
{
  cout << "Cache singleton dtor: AutoSaving" << endl;
  this->AutoSave();

  cout << "Cache singleton dtor: Deleting all cache branches" << endl;
  if(fCacheMap) {
    map<string, CacheBranchI * >::iterator citer;
    for(citer = fCacheMap->begin(); citer != fCacheMap->end(); ++citer) {
      CacheBranchI * branch = citer->second;
      if(branch) {
        delete branch;
        branch = 0;
      }
    }
    fCacheMap->clear();
    delete fCacheMap;
  }
  if(fCacheFile) {
    fCacheFile->Close();
    delete fCacheFile;
  }
  fInstance = 0;
}
//____________________________________________________________________________
Cache * Cache::Instance()
{
  if(fInstance == 0) {
    static Cache::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Cache;

    fInstance->fCacheMap = new map<string, CacheBranchI * >;
    fInstance->AutoLoad();
  }
  return fInstance;
}
//____________________________________________________________________________
CacheBranchI * Cache::FindCacheBranch(string key)
{
  map<string, CacheBranchI *>::const_iterator map_iter = fCacheMap->find(key);

  if (map_iter == fCacheMap->end()) return 0;
  return map_iter->second;
}
//____________________________________________________________________________
void Cache::AddCacheBranch(string key, CacheBranchI * branch)
{
  fCacheMap->insert( map<string, CacheBranchI *>::value_type(key,branch) );
}
//____________________________________________________________________________
string Cache::CacheBranchKey(string k0, string k1, string k2) const
{
  ostringstream key;

  key << k0;
  if(k1.size()>0) key << "/" << k1;
  if(k2.size()>0) key << "/" << k2;

  return key.str();
}
//____________________________________________________________________________
void Cache::AutoLoad(void)
{
  LOG("Cache", pNOTICE) << "AutoLoading Cache";

  this->OpenCacheFile();

  if(!fCacheFile) return;

  TDirectory * cache = (TDirectory *) fCacheFile->Get("cache");
  if(!cache) {
   LOG("Cache", pNOTICE) << "Loaded cache is empty!";
   return;
  }
  cache->cd();

  TList * keys = (TList*) cache->Get("key_list");
  TIter kiter(keys);
  TObjString * keyobj = 0;

  int ib=0;
  while ((keyobj = (TObjString *)kiter.Next())) {

    string key = string(keyobj->GetString().Data());

    ostringstream bname;
    bname << "buffer_" << ib++;

    CacheBranchI * buffer = (CacheBranchI*) cache->Get(bname.str().c_str());
    if(buffer) {
     fCacheMap->insert( map<string, CacheBranchI *>::value_type(key,buffer) );
    }
  }

  LOG("Cache", pNOTICE) << "Cache loaded...";
  LOG("Cache", pNOTICE) << *this;
}
//____________________________________________________________________________
void Cache::AutoSave(void)
{
  if(!fCacheFile) {
    cout << " no cache file is open!" << endl;
    return;
  }

  fCacheFile->cd();

  TDirectory * cache = dynamic_cast<TDirectory *> (fCacheFile->Get("cache"));
  if(!cache) {
    cache = new TDirectory("cache", "GENIE Cache");
    cache->Write("cache");
  }
  cache->cd();

  int ib=0;
  TList * keys = new TList;
  keys->SetOwner(true);

  map<string, CacheBranchI * >::iterator citer;
  for(citer = fCacheMap->begin(); citer != fCacheMap->end(); ++citer) {
    string key = citer->first;
    CacheBranchI * branch = citer->second;
    if(branch) {
      cout << " saving cache buffer: " << key << endl;

      ostringstream bname;
      bname << "buffer_" << ib++;

      keys->Add(new TObjString(key.c_str()));
      branch->Write(bname.str().c_str(), TObject::kOverwrite);
    }
  }
  keys->Write("key_list", TObject::kSingleKey&&TObject::kOverwrite);

  keys->Clear();
  delete keys;
}
//____________________________________________________________________________
void Cache::OpenCacheFile(void)
{
  fCacheFile = 0;

  if(!gSystem->Getenv("GCACHEFILE") ) {
    LOG("Cache", pWARN)
     << "$GCACHEFILE was not set. Cache buffer loading/saving is disabled";
    return;
  }

  string file = gSystem->Getenv("GCACHEFILE");
  LOG("Cache", pNOTICE) << "$GCACHEFILE env.var = " << file;

  fCacheFile = new TFile(file.c_str(),"update");

  if(!fCacheFile->IsOpen()) {
     delete fCacheFile;
     fCacheFile = 0;
     LOG("Cache", pWARN) << "Could not open cache file: " << file;
  }
}
//____________________________________________________________________________
void Cache::Print(ostream & stream) const
{
  stream << "\n [-] GENIE Cache Buffers:";
  stream << "\n  |";
  map<string, CacheBranchI * >::iterator citer;
  for(citer = fCacheMap->begin(); citer != fCacheMap->end(); ++citer) {
    string key = citer->first;
    CacheBranchI * branch = citer->second;
    stream << "\n  |--o  " << key;
    if(!branch) {
      stream << " *** NULL *** ";
    }
  }
  stream << "\n";
}
//___________________________________________________________________________

} // genie namespace


