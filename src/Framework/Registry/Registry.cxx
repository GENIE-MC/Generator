//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 11, 2007 - CA
   Added 'bool ItemIsLocal(RgKey) const', 'void OverrideGlobalDef(RgKey)', 
   and 'void LinkToGlobalDef(RgKey)' to handle the distinction between 
   global/local configuration items. 
   The templated function 'GetValueOrUseDefault' was modifed to mark items 
   with global status & not return them so that an updated default can be 
   cascaded through the entire pool of instantiated algorithms.
 @ Sep 30, 2009 - CA
   Added 'RgType_t ItemType(RgKey) const', 'RgKeyList FindKeys(RgKey) const'

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFolder.h>
#include <TObjString.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Registry/RegistryItemTypeId.h"

using namespace genie;

using std::setw;
using std::setfill;
using std::istream;
using std::cout;
using std::endl;
using std::ostringstream;

//____________________________________________________________________________
namespace genie {
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 template<class T> void SetRegistryItem(Registry * r, RgKey key, T item)
 {
    string itemtype = typeid(item).name();
    LOG("Registry", pINFO)
             << "Set item [" << itemtype << "]: key = "
                                      << key << " --> value = " << item;
    bool lock = r->ItemIsLocked(key); // store, could be true but inhibited
    RegistryItem<T> * reg_item = new RegistryItem<T>(item,lock);
    RgIMapPair config_entry(key, reg_item);
    r->Set(config_entry);
 }
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 template<class T> T GetValueOrUseDefault(
                          Registry * r, RgKey key, T def, bool set_def)
 {
  // Return the requested registry item. If it does not exist return
  // the input default value (in this case, if set_def is true it can 
  // override a lock and add the input default as a new registry item)

   T value;
   if(r->Exists(key)) { 
      if(r->ItemIsLocal(key)) {
         r->Get(key,value); 
         return value;
      }
   }
   value = def;
   bool was_locked = r->IsLocked();
   if(was_locked) r->UnLock();

   if(set_def) {
      r->Set(key, value);
      r->LinkToGlobalDef(key);
   }
   if(was_locked) r->Lock();
   return value;
 }
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ostream & operator << (ostream & stream, const Registry & registry)
 {
   registry.Print(stream);
   return stream;
 }
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//____________________________________________________________________________
Registry::Registry() 
{
  this->Init();
}
//____________________________________________________________________________
Registry::Registry(string name, bool isReadOnly) :
fName             ( name       ),
fIsReadOnly       ( isReadOnly ),
fInhibitItemLocks ( false      )
{

}
//____________________________________________________________________________
Registry::Registry(const Registry & registry) :
fName("uninitialised"),
fIsReadOnly(false)
{
  this->Copy(registry);
}
//____________________________________________________________________________
Registry::~Registry()
{
  this->Clear(true);
}
//____________________________________________________________________________
void Registry::operator() (RgKey key, RgInt item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (RgKey key, RgBool item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (RgKey key, RgDbl item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (RgKey key, RgCChAr item)
{
  RgStr item2 = RgStr(item); // "const char *" -> "string"
  this->Set(key,item2);
}
//____________________________________________________________________________
void Registry::operator() (RgKey key, RgStr item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
Registry & Registry::operator = (const Registry & reg)
{
  this->Copy(reg);
  return (*this);
}
//____________________________________________________________________________
Registry & Registry::operator += (const Registry & reg)
{
  this->Append(reg);
  return (*this);
}
//____________________________________________________________________________
void Registry::Lock(void)
{
  fIsReadOnly = true;
}
//____________________________________________________________________________
void Registry::UnLock(void)
{
  fIsReadOnly = false;
}
//____________________________________________________________________________
bool Registry::IsLocked(void) const
{
  return fIsReadOnly;
}
//____________________________________________________________________________
void Registry::InhibitItemLocks(void)
{
  fInhibitItemLocks = true;
}
//____________________________________________________________________________
void Registry::RestoreItemLocks(void)
{
  fInhibitItemLocks = false;
}
//____________________________________________________________________________
bool Registry::ItemLocksAreActive(void) const
{
  return !fInhibitItemLocks;
}
//____________________________________________________________________________
bool Registry::ItemIsLocal(RgKey key) const
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     bool is_local = entry->second->IsLocal();
     return is_local;
  } else {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("Registry", pDEBUG)
        << "*** Was asked to check 'local' flag on non-existing item: [" 
        << key << "]";
#endif
  }
  return false;
}
//____________________________________________________________________________
void Registry::OverrideGlobalDef(RgKey key)
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     entry->second->SetLocal(true);
  } else {
     LOG("Registry", pWARN)
        << "*** Can't give 'local' status to  non-existem item [" 
        << key << "]"; 
  }
}
//____________________________________________________________________________
void Registry::LinkToGlobalDef(RgKey key)
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     entry->second->SetLocal(false);
  } else {
     LOG("Registry", pWARN)
        << "*** Can't give 'global' status to  non-existem item [" 
        << key << "]"; 
  }
}
//____________________________________________________________________________
bool Registry::ItemIsLocked(RgKey key) const
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     bool is_locked = entry->second->IsLocked();
     return is_locked;
  } else {
/*
     LOG("Registry", pDEBUG)
        << "*** Was asked to check lock on non-existing item: [" 
        << key << "]";
*/
  }  
  return false;
}
//____________________________________________________________________________
void Registry::LockItem(RgKey key)
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     entry->second->Lock();
  } else {
     LOG("Registry", pWARN)
           << "*** Can't lock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
void Registry::UnLockItem(RgKey key)
{
  if( this->Exists(key) ) {
     RgIMapConstIter entry = fRegistry.find(key);
     entry->second->UnLock();
  } else {
    LOG("Registry", pWARN)
       << "*** Can't unlock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
bool Registry::CanSetItem(RgKey key) const
{
  bool locked_item       = this->ItemIsLocked(key);
  bool active_item_locks = this->ItemLocksAreActive();
  bool locked_registry   = this->IsLocked();

  bool can_set = !locked_registry && ( !locked_item || !active_item_locks );

  return can_set;
}
//____________________________________________________________________________
void Registry::Set(RgIMapPair entry)
{
  RgKey key = entry.first;
  if( this->CanSetItem(key) ) {
    this->DeleteEntry(key);
    fRegistry.insert(entry);
  } else {
     LOG("Registry", pWARN)
             << "*** Registry item [" << key << "] can not be set";
  }
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgBool item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgInt item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgDbl item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgCChAr item)
{
  RgStr item2 = RgStr(item); // "const char *" -> "string"
  this->Set(key, item2);
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgStr item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgAlg item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgH1F item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgH2F item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(RgKey key, RgTree item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Get(RgKey key, const RegistryItemI * & item) const
{
  RgIMapConstIter entry = this->SafeFind(key);
  item = entry->second;
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgBool & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Get an RgBool item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgBool> * ri = dynamic_cast<RegistryItem<RgBool>*> (rib);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
#endif
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgInt & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgInt item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgInt> * ri = dynamic_cast< RegistryItem<RgInt> * > (rib);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
#endif
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgDbl & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgDbl item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgDbl> * ri = dynamic_cast<RegistryItem<RgDbl>*> (rib);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
#endif
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgStr & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgStr item with  key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgStr> * ri = dynamic_cast<RegistryItem<RgStr>*> (rib);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
#endif
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgAlg & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgAlg item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgAlg> * ri = dynamic_cast<RegistryItem<RgAlg>*> (rib);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
#endif

  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgH1F & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgH1F item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgH1F> *ri = dynamic_cast<RegistryItem<RgH1F>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TH1F param = " << key;
  }
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgH2F & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgH2F item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgH2F> *ri = dynamic_cast<RegistryItem<RgH2F>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TH2F param = " << key;
  }
}
//____________________________________________________________________________
void Registry::Get(RgKey key, RgTree & item) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Registry", pDEBUG) << "Getting an RgTree item with key: " << key;
#endif

  RgIMapConstIter entry = this->SafeFind(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgTree> *ri =  dynamic_cast<RegistryItem<RgTree>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TTree param = " << key;
  }
}
//____________________________________________________________________________
RgBool Registry::GetBool(RgKey key) const
{
  RgBool value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
RgInt Registry::GetInt(RgKey key) const
{
  RgInt value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
RgDbl Registry::GetDouble(RgKey key) const
{
  RgDbl value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
RgStr Registry::GetString(RgKey key) const
{
  RgStr value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
RgAlg Registry::GetAlg(RgKey key) const
{
  RgAlg value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
RgH1F Registry::GetH1F(RgKey key) const
{
  RgIMapConstIter entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgH1F> *ri = dynamic_cast<RegistryItem<RgH1F>*> (rib);

  RgH1F item = ri->Data();
  return item;
}
//____________________________________________________________________________
RgH2F Registry::GetH2F(RgKey key) const
{
  RgIMapConstIter entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgH2F> *ri = dynamic_cast<RegistryItem<RgH2F>*> (rib);

  RgH2F item = ri->Data();
  return item;
}
//____________________________________________________________________________
RgTree Registry::GetTree(RgKey key) const
{
  RgIMapConstIter entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<RgTree> *ri = dynamic_cast<RegistryItem<RgTree>*> (rib);

  RgTree item = ri->Data();
  return item;
}
//____________________________________________________________________________
RgBool Registry::GetBoolDef(RgKey key, RgBool def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
int Registry::GetIntDef(RgKey key, int def_opt, bool set_def)
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
double Registry::GetDoubleDef(RgKey key, double def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
string Registry::GetStringDef(RgKey key, string def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
RgAlg Registry::GetAlgDef(RgKey key, RgAlg  def_opt, bool set_def)
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
RgIMapConstIter Registry::SafeFind(RgKey key) const
{
  RgIMapConstIter entry = fRegistry.find(key);
  if (entry!=fRegistry.end()) {
    return entry;
  }
  LOG("Registry/SafeFind", pFATAL)
       << "*** Key: " << key
       << " does not exist in registry: " << this->Name();
  gAbortingInErr = true;
  exit(1);    
}
//____________________________________________________________________________
bool Registry::Exists(RgKey key) const
{
  RgIMapConstIter entry = fRegistry.find(key);
  return (entry!=fRegistry.end());
}
//____________________________________________________________________________
bool Registry::DeleteEntry(RgKey key)
{
  if(!fIsReadOnly && Exists(key)) {
      RgIMapIter entry = fRegistry.find(key);
      RegistryItemI * item = entry->second;
      delete item;
      item = 0;
      fRegistry.erase(entry);
      return true;
  }
  return false;
}
//____________________________________________________________________________
int Registry::NEntries(void) const
{
  RgIMapSizeType reg_size = fRegistry.size();
  return (const int) reg_size;
}
//____________________________________________________________________________
void Registry::SetName(string name)
{
  if(! fIsReadOnly) fName = name;
  else {
     LOG("Registry", pWARN)
        << "*** Registry is locked - Can not change its name";
  }
}
//____________________________________________________________________________
string Registry::Name(void) const
{
  return fName;
}
//____________________________________________________________________________
void Registry::AssertExistence(RgKey key0) const
{
  if ( ! this->Exists(key0) ) {
     LOG("Registry", pERROR) << (*this);
     LOG("Registry", pFATAL)
           << "*** Key: " << key0
             << " does not exist in registry: " << this->Name();
     exit(1);    
  }
}
//____________________________________________________________________________
void Registry::AssertExistence(RgKey key0, RgKey key1) const
{
  this->AssertExistence(key0);
  this->AssertExistence(key1);
}
//____________________________________________________________________________
void Registry::AssertExistence(RgKey key0, RgKey key1, RgKey key2) const
{
  this->AssertExistence(key0);
  this->AssertExistence(key1);
  this->AssertExistence(key2);
}
//____________________________________________________________________________
void Registry::CopyToFolder(TFolder * folder) const
{
  LOG("Registry", pINFO) << "Converting Registry to TFolder";

  folder->SetOwner(true);

  RgIMapConstIter reg_iter;

  for(reg_iter = this->fRegistry.begin();
                          reg_iter != this->fRegistry.end(); reg_iter++) {

     ostringstream   entry;
     string          key   = reg_iter->first;
     RegistryItemI * ritem = reg_iter->second;

     RgType_t type  = ritem->TypeInfo();
     string   stype = RgType::AsString(type);

     entry << "key:" << key << ";type:" << stype;

     if(type == kRgBool) {
        entry << ";value: " << this->GetBool(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == kRgDbl) {
        entry << ";value: " << this->GetDouble(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == kRgInt) {
        entry << ";value: " << this->GetInt(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == kRgStr) {
        entry << ";value: " << this->GetString(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == kRgAlg) {
        entry << ";value: " << this->GetAlg(key).name 
              << "/" << this->GetAlg(key).config;
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));

     } else if (type == kRgH1F) {
     } else if (type == kRgH2F) {
     } else if (type == kRgTree) {
     } else {}
  }// registry iterator
}
//____________________________________________________________________________
void Registry::Print(ostream & stream) const
{
// Prints the registry to the specified stream
//
   stream << endl;
   stream << "[-] Registry name: [" << Name() << "]";

   stream << " - Write Status: ";
   if(fIsReadOnly)       { stream << "[locked]";   }
   else                  { stream << "[unlocked]"; }

   stream << " - Inhibited Item Locking: ";
   if(fInhibitItemLocks) { stream << "[on]";       }
   else                  { stream << "[off]";      }

   stream << " - # entries: " << setfill(' ') << setw(3) << fRegistry.size()
          << endl;

   RgIMapConstIter rcit = fRegistry.begin();
   for( ; rcit != fRegistry.end(); rcit++) {

     RgKey           key   = rcit->first;
     RegistryItemI * ritem = rcit->second;
     if(ritem) {
        RgType_t type  = ritem->TypeInfo();
        string   stype = RgType::AsString(type);

        string key_lbl  = string("> ") + key;
        string type_lbl = string("[") + stype + string("] ");
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("Registry", pDEBUG)
                      << "Printing [" << stype << "] item named = " << key;
#endif
        stream << " |" << setfill('-') << setw(50) << key_lbl
               << setfill(' ') << setw(10) << type_lbl;
        ritem->Print(stream);
        stream << endl;
     } else {
        LOG("Registry", pERROR) << "Null RegistryItemI with key = " << key;
     }
   }// registry iterator
}
//____________________________________________________________________________
void Registry::Copy(const Registry & registry)
{
// Copies the input registry
//
  LOG("Registry", pINFO) 
         << "Copying registry " << registry.Name() << " to " << this->Name();

  if(this->IsLocked()) {
   LOG("Registry", pWARN) << "Registry is locked. Can't copy input entries!";
   return;
  }

  this->Init();
  this->Clear();
  this->Append(registry);

  fName             = registry.fName;
  fIsReadOnly       = registry.fIsReadOnly;
  fInhibitItemLocks = registry.fInhibitItemLocks;
}
//____________________________________________________________________________
void Registry::Append(const Registry & registry, RgKey prefix)
{
// Appends the input registry entries (& their locks)

  LOG("Registry", pINFO) 
       << "Appending registry " << registry.Name() << " to " << this->Name();

  if(this->IsLocked()) {
   LOG("Registry", pWARN) << "Registry is locked. Can't copy input entries!";
   return;
  }

  this->InhibitItemLocks();

  RgIMapConstIter reg_iter;
  for(reg_iter = registry.fRegistry.begin();
                      reg_iter != registry.fRegistry.end(); reg_iter++) {

     RgKey name     = reg_iter->first;
     RgKey new_name = prefix + name;

     if ( fRegistry.count( new_name ) > 0 ) continue ;

     RgType_t type  = reg_iter -> second -> TypeInfo();
     string   stype = RgType::AsString(type);

     LOG("Registry", pINFO)
         << "Copying [" << stype << "] item named = " 
                                         << name << " as " << new_name;

     RegistryItemI * cri = registry.CloneRegistryItem( name ) ; // cloned registry item

     RgIMapPair reg_entry(new_name, cri);

     if ( ! fRegistry.insert(reg_entry).second ) {
    	 // The registry already contained an entry with key new_name
    	 //   so the new registryItem has to be deleted or we leak memory.
    	 //   This should not happened as a check is performed
    	 LOG("Registry", pERROR ) << "Failing to insert item " << new_name ;
    	 delete cri ;
     }
   } // loop on the incoming registry items
}
//____________________________________________________________________________
void Registry::Merge(const Registry & registry, RgKey prefix)
{
// Add the input registry entries (& their locks)
// and updated entries already present

  LOG("Registry", pINFO)
       << "Appending registry " << registry.Name() << " to " << this->Name();

  if(this->IsLocked()) {
   LOG("Registry", pWARN) << "Registry is locked. Can't copy input entries!";
   return;
  }

  this->InhibitItemLocks();

  RgIMapConstIter reg_iter;
  for(reg_iter = registry.fRegistry.begin();
                      reg_iter != registry.fRegistry.end(); reg_iter++) {

     RgKey name     = reg_iter->first;
     RgKey new_name = prefix + name;

     RgType_t type  = reg_iter -> second -> TypeInfo();
     string   stype = RgType::AsString(type);

     LOG("Registry", pINFO)
         << "Copying [" << stype << "] item named = "
                                         << name << " as " << new_name;

     RegistryItemI * cri = registry.CloneRegistryItem( name ) ; // cloned registry item

     if ( fRegistry.count( new_name ) > 0 ) {

       RegistryItemI * old_ri = fRegistry[new_name] ;
   	   delete old_ri ;
     }

     fRegistry[new_name] = cri ;

   } // loop on the incoming registry items

}//____________________________________________________________________________
RgType_t Registry::ItemType(RgKey key) const
{
  RgIMapConstIter reg_iter = fRegistry.find(key);
  if(reg_iter != fRegistry.end()) {
     RegistryItemI * ri = reg_iter->second;
     RgType_t type  = ri->TypeInfo();
     return type;
  }
  return kRgUndefined;
}
//____________________________________________________________________________
RgKeyList Registry::FindKeys(RgKey key_part) const
{
  RgKeyList klist;

  RgIMapConstIter reg_iter = fRegistry.begin();
  for( ; reg_iter != fRegistry.end(); reg_iter++) {
    RgKey key = reg_iter->first;
    if (key.find(key_part) != string::npos) {
      klist.push_back(key);
    }
  }

  return klist;
}
//____________________________________________________________________________
void Registry::Init(void)
{
// initialize registry properties

  fName              = "NoName";
  fIsReadOnly        = false;
  fInhibitItemLocks  = false;
}
//____________________________________________________________________________
void Registry::Clear(bool force)
{
// clean all registry entries

  LOG("Registry", pINFO) 
      << "Cleaning-up [force unlock = " << ((force)?"true":"false")
      << "] registry: " << this->Name();
  if(!force) {
   if(this->IsLocked()) {
      LOG("Registry", pWARN) << "Registry is locked. Can't clear its entries";
      return;
    }
  }
  RgIMapIter rit;
  for(rit = fRegistry.begin(); rit != fRegistry.end(); rit++) {
     RgKey           name = rit->first;
     RegistryItemI * item = rit->second;
     if(!item) {
       LOG("Registry", pWARN) << "Item with key = " << name << " is null!";
     }
     delete item;
     item = 0;
  }
  fRegistry.clear();
}
//____________________________________________________________________________
RegistryItemI * Registry::CloneRegistryItem( const RgKey & key ) const {

	std::map<RgKey, RegistryItemI*>::const_iterator it = fRegistry.find( key ) ;

	if ( it == fRegistry.end() ) {
		LOG("Registry", pFATAL) << "Item " << key << " not found while cloning for registry " << Name() ;
		exit( 0 ) ;
	}

	RegistryItemI * ri  = it -> second ;

     bool     ilk   = ri->IsLocked();
     RgType_t type  = ri->TypeInfo();
     string   stype = RgType::AsString(type);

     RegistryItemI * cri = 0; // cloned registry item
     if (type == kRgBool)
           cri = new RegistryItem<RgBool>( GetBool(key), ilk);
     else if (type == kRgDbl)
           cri = new RegistryItem<RgDbl> ( GetDouble(key), ilk);
     else if (type == kRgInt)
           cri = new RegistryItem<RgInt> ( GetInt(key), ilk);
     else if (type == kRgStr)
           cri = new RegistryItem<RgStr> ( GetString(key), ilk);
     else if (type == kRgAlg)
           cri = new RegistryItem<RgAlg> ( GetAlg(key), ilk);
     else if (type == kRgH1F) {
           RgH1F histo = GetH1F(key);
           if(histo) {
               RgH1F chisto = new TH1F(*histo);
               LOG("Registry", pDEBUG) << chisto->GetName();
               cri = new RegistryItem<RgH1F>(chisto,ilk);
           } else {
             LOG("Registry", pERROR)
               << "Null TH1F with key = " << key << " - not copied";
           }
     } else if (type == kRgH2F) {
           RgH2F histo = GetH2F(key);
           if(histo) {
               RgH2F chisto = new TH2F(*histo);
               LOG("Registry", pDEBUG) << chisto->GetName();
               cri = new RegistryItem<RgH2F>(chisto,ilk);
           } else {
             LOG("Registry", pERROR)
               << "Null TH2F with key = " << key << " - not copied";
           }
     } else if (type == kRgTree) {
           RgTree tree = GetTree(key);
           if(tree) {
               //TTree * ctree = new TTree(*tree);
               TTree * ctree = tree->CopyTree("1");
               LOG("Registry", pDEBUG) << ctree->GetName();
               cri = new RegistryItem<RgTree>(ctree,ilk);
           } else {
             LOG("Registry", pERROR)
               << "Null TTree with key = " << key << " - not copied";
           }
     } else {

       LOG( "Registry", pFATAL ) << "Item " << key << " not cloned because its type is not implemented " ;
       exit( 0 ) ;
     }

     return cri ;

}

