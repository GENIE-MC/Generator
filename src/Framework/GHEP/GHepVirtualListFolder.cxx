//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - July 16, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepVirtualListFolder.h"
#include "Framework/GHEP/GHepVirtualList.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
GHepVirtualListFolder * GHepVirtualListFolder::fInstance = 0;
//____________________________________________________________________________
GHepVirtualListFolder::GHepVirtualListFolder()
{
  fInstance =  0;
}
//____________________________________________________________________________
GHepVirtualListFolder::~GHepVirtualListFolder()
{
  this->Clear();
  fInstance = 0;
}
//____________________________________________________________________________
GHepVirtualListFolder * GHepVirtualListFolder::Instance()
{
  if(fInstance == 0) {

    static GHepVirtualListFolder::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new GHepVirtualListFolder;
  }
  return fInstance;
}
//____________________________________________________________________________
void GHepVirtualListFolder::AddToVirtualList(string listname, GHepParticle* p)
{
// Adds a particle to the named virtual list - if the list does not exists it
// creates it first.
// The virtual list has no ownership of its
  bool exists = this->VirtualListExists(listname);

  if(!exists) this->AddVirtualList(listname);

  int n = fVirtualListMap[listname]->GetEntries();
  GHepVirtualList * vl = fVirtualListMap[listname];
  (*vl)[n] = (TObject*)p;
}
//____________________________________________________________________________
bool GHepVirtualListFolder::VirtualListExists(string listname)
{
// checks whether a virtual list exists

  bool exists = (fVirtualListMap.count(listname) == 1);
  return exists;
}
//____________________________________________________________________________
void GHepVirtualListFolder::RemoveList(string listname)
{
// removes the input virtual list (if it exists)

  bool exists = (fVirtualListMap.count(listname) == 1);
  if(!exists) return;

  map<string, GHepVirtualList *>::iterator
                                vlmiter = fVirtualListMap.find(listname);

  GHepVirtualList * vlist = vlmiter->second;
  if(vlist) delete vlist;
  vlist = 0;

  fVirtualListMap.erase(listname);
}
//____________________________________________________________________________
void GHepVirtualListFolder::Clear(void)
{
// removes all virtual lists

  map<string, GHepVirtualList *>::iterator vlmiter;
  for (vlmiter = fVirtualListMap.begin();
                          vlmiter != fVirtualListMap.end(); ++vlmiter) {
     GHepVirtualList * vlist = vlmiter->second;
     if(vlist) delete vlist;
     vlist = 0;
  }
  fVirtualListMap.clear();
}
//____________________________________________________________________________
GHepVirtualList * GHepVirtualListFolder::VirtualList(string listname)
{
  bool exists = this->VirtualListExists(listname);

  if(!exists) return 0;
  else        return fVirtualListMap[listname];
}
//____________________________________________________________________________
void GHepVirtualListFolder::AddVirtualList(string listname)
{
  GHepVirtualList * vl = new GHepVirtualList;
  fVirtualListMap.insert(
                   map<string, GHepVirtualList *>::value_type(listname,vl) );
}
//____________________________________________________________________________
