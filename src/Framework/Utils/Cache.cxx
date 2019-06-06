//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 24, 2013 - CA
   Cache is not autoloaded and use of variables $GCACHEFILE is no longer
   supported. Instead, call Cache::OpenCacheFile(string filename) explicitly.
   Now cached data are stored in the top-level 'directory'.
*/
//____________________________________________________________________________

#include <sstream>
#include <iostream>

#include <TSystem.h>
#include <TDirectory.h>
#include <TList.h>
#include <TObjString.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchI.h"

using std::ostringstream;
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
  this->Save();

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
void Cache::RmCacheBranch(string key)
{
  LOG("Cache", pNOTICE) << "Removing cache branch: " << key;

}
//____________________________________________________________________________
void Cache::RmAllCacheBranches(void)
{
  LOG("Cache", pNOTICE) << "Removing cache branches";

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
  }
}
//____________________________________________________________________________
void Cache::RmMatchedCacheBranches(string key_substring)
{
  LOG("Cache", pNOTICE) << "Removing cache branches: *"<< key_substring<< "*";

}
//____________________________________________________________________________
void Cache::Load(void)
{
  LOG("Cache", pNOTICE) << "Loading cache";

  if(!fCacheFile) return;
  TList * keys = (TList*) fCacheFile->Get("key_list");
  TIter kiter(keys);
  TObjString * keyobj = 0;
  int ib=0;
  while ((keyobj = (TObjString *)kiter.Next())) {
    string key = string(keyobj->GetString().Data());
    ostringstream bname;
    bname << "buffer_" << ib++;
    CacheBranchI * buffer = (CacheBranchI*) fCacheFile->Get(bname.str().c_str());
    if(buffer) {
     fCacheMap->insert( map<string, CacheBranchI *>::value_type(key,buffer) );
    }
  }
  LOG("Cache", pNOTICE) << "Cache loaded...";
  LOG("Cache", pNOTICE) << *this;
}
//____________________________________________________________________________
void Cache::Save(void)
{
  if(!fCacheFile) {
    return;
  }
  fCacheFile->cd();

  int ib=0;
  TList * keys = new TList;
  keys->SetOwner(true);

  map<string, CacheBranchI * >::iterator citer;
  for(citer = fCacheMap->begin(); citer != fCacheMap->end(); ++citer) {
    string key = citer->first;
    CacheBranchI * branch = citer->second;
    if(branch) {
      ostringstream bname;
      bname << "buffer_" << ib++;
      keys->Add(new TObjString(key.c_str()));
      branch->Write(bname.str().c_str(), TObject::kOverwrite);
    }
  }
  keys->Write("key_list", TObject::kSingleKey | TObject::kOverwrite );

  keys->Clear();
  delete keys;
}
//____________________________________________________________________________
void Cache::OpenCacheFile(string filename)
{
  if(filename.size() == 0) return;

  if(fCacheFile) {
    if(fCacheFile->IsOpen()) {
       delete fCacheFile;
       fCacheFile = 0;
    }
  }

  LOG("Cache", pNOTICE) << "Using cache file: " << filename;

  fCacheFile = new TFile(filename.c_str(),"update");
  if(!fCacheFile->IsOpen()) {
     delete fCacheFile;
     fCacheFile = 0;
     LOG("Cache", pWARN) << "Could not open cache file: " << filename;
  }

  this->Load();
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


