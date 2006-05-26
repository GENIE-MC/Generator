//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <iomanip>

#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFolder.h>
#include <TObjString.h>

#include "Messenger/Messenger.h"
#include "Registry/Registry.h"

using namespace genie;

using std::setw;
using std::setfill;
using std::istream;
using std::endl;
using std::ostringstream;

//____________________________________________________________________________
namespace genie {

 template<class T> void SetRegistryItem(Registry * r, string key, T item)
 {
    string itemtype = typeid(item).name();
    LOG("Registry", pDEBUG)
             << "Set item [" << itemtype << "]: key = "
                                      << key << " --> value = " << item;
    bool lock = r->ItemIsLocked(key); // store, could be true but inhibited
    RegistryItem<T> * reg_item = new RegistryItem<T>(item,lock);
    pair<string, RegistryItemI *> config_entry(key, reg_item);
    r->Set(config_entry);
 }
 //...........................................................................
 template<class T> T GetValueOrUseDefault(
                          Registry * r, string key, T def, bool set_def)
 {
  // Return the requested registry item. If it does not exist return
  // the input default value (in this case, if set_def is true it can 
  // override a lock and add the input default as a new registry item)

   T value;
   if(r->Exists(key)) { 
      r->Get(key,value); return value;
   }

   value = def;
   bool was_locked = r->IsLocked();
   if(was_locked) r->UnLock();

   if(set_def) r->Set(key, value);
   if(was_locked) r->Lock();

   return value;
 }
 //...........................................................................
 ostream & operator << (ostream & stream, const Registry & registry)
 {
   registry.Print(stream);
   return stream;
 }
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
void Registry::operator() (string key, int item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (string key, bool item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (string key, double item)
{
  this->Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (string key, const char * item)
{
  this->Set(key, string(item));
}
//____________________________________________________________________________
void Registry::operator() (string key, string item)
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
  if( this->IsLocked() ) {
      LOG("Registry", pWARN)
               << "Registry is locked. Can not add more entries!";
      return (*this);
  }

  Registry inp_reg(reg);

  map<string, RegistryItemI *>::const_iterator reg_iter;
  for(reg_iter = inp_reg.fRegistry.begin();
                   reg_iter != inp_reg.fRegistry.end(); reg_iter++) {
     string          name = reg_iter->first;
     RegistryItemI * item = reg_iter->second;

     pair<string, RegistryItemI *> reg_entry(name, item);
     this->fRegistry.insert(reg_entry);
  }
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
bool Registry::ItemIsLocked(string key) const
{
  if( this->Exists(key) ) {

     map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);

     bool is_locked = entry->second->IsLocked();
     return is_locked;
  }
  return false;
}
//____________________________________________________________________________
void Registry::LockItem(string key)
{
  if( this->Exists(key) ) {
     map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
     entry->second->Lock();
  } else {
     LOG("Registry", pWARN)
           << "*** Can't lock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
void Registry::UnLockItem(string key)
{
  if( this->Exists(key) ) {
     map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
     entry->second->UnLock();
  } else {
    LOG("Registry", pWARN)
       << "*** Can't unlock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
bool Registry::CanSetItem(string key) const
{
  bool locked_item       = this->ItemIsLocked(key);
  bool active_item_locks = this->ItemLocksAreActive();
  bool locked_registry   = this->IsLocked();

  bool can_set = !locked_registry && ( !locked_item || !active_item_locks );

  return can_set;
}
//____________________________________________________________________________
void Registry::Set(pair<string, RegistryItemI *> entry)
{
  string key = entry.first;

  if( this->CanSetItem(key) ) {
    this->DeleteEntry(key);
    fRegistry.insert(entry);
  } else {
     LOG("Registry", pWARN)
             << "*** Registry item [" << key << "] can not be set";
  }
}
//____________________________________________________________________________
void Registry::Set(string key, bool item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, int item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, double item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, const char* item)
{
  this->Set(key, string(item));
}
//____________________________________________________________________________
void Registry::Set(string key, string item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, TH1F * item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, TH2F * item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Set(string key, TTree * item)
{
  SetRegistryItem(this, key, item); // call templated set method
}
//____________________________________________________________________________
void Registry::Get(string key, const RegistryItemI * item) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  item = entry->second;
}
//____________________________________________________________________________
void Registry::Get(string key, bool & item) const
{
  LOG("Registry", pDEBUG) << "Get item [bool]: key = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<bool> * ri = dynamic_cast<RegistryItem<bool>*> (rib);

  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(string key, int & item) const
{
  LOG("Registry", pDEBUG) << "Getting item (int) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<int> * ri = dynamic_cast< RegistryItem<int> * > (rib);

  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(string key, double & item) const
{
  LOG("Registry", pDEBUG) << "Getting item (double) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<double> * ri = dynamic_cast<RegistryItem<double>*> (rib);

  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(string key, string & item) const
{
  LOG("Registry", pDEBUG) << "Getting item (string) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<string> * ri = dynamic_cast<RegistryItem<string>*> (rib);

  LOG("Registry", pDEBUG) << "Item value = " << ri->Data();
  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(string key, TH1F * item) const
{
  LOG("Registry", pDEBUG) << "Getting item (TH1F) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TH1F*> *ri = dynamic_cast<RegistryItem<TH1F*>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TH1F param = " << key;
  }
}
//____________________________________________________________________________
void Registry::Get(string key, TH2F * item) const
{
  LOG("Registry", pDEBUG) << "Getting item (TH2F) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TH2F*> *ri = dynamic_cast<RegistryItem<TH2F*>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TH2F param = " << key;
  }
}
//____________________________________________________________________________
void Registry::Get(string key, TTree * item) const
{
  LOG("Registry", pDEBUG) << "Getting item (TTree) = " << key;

  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TTree*> *ri =  dynamic_cast<RegistryItem<TTree*>*> (rib);
  item = ri->Data();

  if(!item) {
    LOG("Registry", pWARN) << "Returned NULL ptr for TTree param = " << key;
  }
}
//____________________________________________________________________________
bool Registry::GetBool(string key) const
{
  bool value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
int Registry::GetInt(string key) const
{
  int value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
double Registry::GetDouble(string key) const
{
  double value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
string Registry::GetString(string key) const
{
  string value;
  this->Get(key, value);
  return value;
}
//____________________________________________________________________________
TH1F * Registry::GetTH1F(string key) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TH1F*> *ri = dynamic_cast<RegistryItem<TH1F*>*> (rib);

  TH1F* item = ri->Data();
  return item;
}
//____________________________________________________________________________
TH2F * Registry::GetTH2F(string key) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TH2F*> *ri = dynamic_cast<RegistryItem<TH2F*>*> (rib);

  TH2F* item = ri->Data();
  return item;
}
//____________________________________________________________________________
TTree * Registry::GetTTree(string key) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
  RegistryItemI * rib = entry->second;
  RegistryItem<TTree*> *ri = dynamic_cast<RegistryItem<TTree*>*> (rib);

  TTree* item = ri->Data();
  return item;
}
//____________________________________________________________________________
bool Registry::GetBoolDef(string key, bool def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
int Registry::GetIntDef(string key, int def_opt, bool set_def)
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
double Registry::GetDoubleDef(string key, double def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
string Registry::GetStringDef(string key, string def_opt, bool set_def) 
{
  return GetValueOrUseDefault(this, key, def_opt, set_def);
}
//____________________________________________________________________________
bool Registry::Exists(string key) const
{
  if (fRegistry.count(key) == 1) return true;
  else                           return false;
}
//____________________________________________________________________________
bool Registry::DeleteEntry(string key)
{
  if(!fIsReadOnly && Exists(key)) {
      map<string, RegistryItemI *>::iterator entry = fRegistry.find(key);
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
  map<string, RegistryItemI *>::size_type reg_size = fRegistry.size();

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
void Registry::AssertExistence(string key0) const
{
  if ( ! this->Exists(key0) ) {
     LOG("Registry", pERROR)
           << "*** Key: " << key0
             << " does not exist in registry: " << this->Name();
  }
  assert( this->Exists(key0) );
}
//____________________________________________________________________________
void Registry::AssertExistence(string key0, string key1) const
{
  this->AssertExistence(key0);
  this->AssertExistence(key1);
}
//____________________________________________________________________________
void Registry::AssertExistence(string key0, string key1, string key2) const
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

  map<string, RegistryItemI *>::const_iterator reg_iter;

  for(reg_iter = this->fRegistry.begin();
                          reg_iter != this->fRegistry.end(); reg_iter++) {

     ostringstream   entry;
     string          key   = reg_iter->first;
     RegistryItemI * ritem = reg_iter->second;
     string          type  = string(ritem->TypeInfo().name());

     entry << "key:" << key << ";type:" << type;

     if(type == "b") {
        entry << ";value: " << this->GetBool(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == "d") {
        entry << ";value: " << this->GetDouble(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == "i") {
        entry << ";value: " << this->GetInt(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == "Ss") {
        entry << ";value: " << this->GetString(key);
        LOG("Registry", pINFO) << "entry = " << entry.str();
        folder->Add(new TObjString(entry.str().c_str()));
     }
     else if (type == "P4TH1F")
     {
     } else if (type == "P4TH2F") {
     } else if (type == "P4TTree") {
     } else {}
  }// registry iterator
}
//____________________________________________________________________________
void Registry::Print(ostream & stream) const
{
   stream << endl;
   stream << "[-] REGISTRY: [" << Name() << "]";

   stream << " - Write Status: ";
   if(fIsReadOnly)       { stream << "[locked]";   }
   else                  { stream << "[unlocked]"; }

   stream << " - Inhibited Item Locking: ";
   if(fInhibitItemLocks) { stream << "[on]";       }
   else                  { stream << "[off]";      }

   stream << " - # entries: " << setfill(' ') << setw(3) << fRegistry.size();

   stream << endl;

   map<string, RegistryItemI *>::const_iterator r_iter;

   for(r_iter = fRegistry.begin(); r_iter != fRegistry.end(); r_iter++) {

     string key = r_iter->first;
     RegistryItemI * ritem = r_iter->second;

     if(ritem) {
        string var_name = string("> ") + key;
        string var_type = string("[") + string(ritem->TypeInfo().name())
                        + string("] ");
        LOG("Registry", pDEBUG)
                    << "Printing [" << var_type << "] item named = " << key;
        stream << " |" << setfill('-') << setw(45) << var_name
               << setfill('.') << setw(12) << var_type;
        ritem->Print(stream);

     } else {
        LOG("Registry", pERROR) << "Null RegistryItemI with key = " << key;
     }
   }// registry iterator
}
//____________________________________________________________________________
void Registry::Copy(const Registry & registry)
{
// copy input registry

  if(this->IsLocked()) {
   LOG("Registry", pWARN) << "Registry is locked. Can't copy input entries!";
   return;
  }

  this->Init();
  this->Clear();

  this->InhibitItemLocks();

  map<string, RegistryItemI *>::const_iterator reg_iter;
  for(reg_iter = registry.fRegistry.begin();
                      reg_iter != registry.fRegistry.end(); reg_iter++) {

     string         name = reg_iter->first;
     RegistryItemI * ri  = reg_iter->second;

     string var_type  = ri->TypeInfo().name();
     bool   item_lock = ri->IsLocked();

     LOG("Registry", pDEBUG)
                << "Copying [" << var_type << "] item named = " << name;

     RegistryItemI * cri = 0; // cloned registry item
     if (var_type == "b" )
           cri = new RegistryItem<bool>
                                (registry.GetBool(name),item_lock);
     else if (var_type == "d" )
           cri = new RegistryItem<double>
                                (registry.GetDouble(name),item_lock);
     else if (var_type == "i" )
           cri = new RegistryItem<int>
                                (registry.GetInt(name),item_lock);
     else if (var_type == "Ss")
           cri = new RegistryItem<string>
                                (registry.GetString(name),item_lock);
     else if (var_type == "P4TH1F") {
           TH1F * histo = registry.GetTH1F(name);
           if(histo) {
               TH1F * chisto = new TH1F(*histo);
               LOG("Registry", pDEBUG) << chisto->GetName();
               cri = new RegistryItem<TH1F*>(chisto,item_lock);
           } else {
             LOG("Registry", pERROR)
               << "Null TH1F with key = " << name << " - not copied";
           }
     } else if (var_type == "P4TH2F") {
           TH2F * histo = registry.GetTH2F(name);
           if(histo) {
               TH2F * chisto = new TH2F(*histo);
               LOG("Registry", pDEBUG) << chisto->GetName();
               cri = new RegistryItem<TH2F*>(chisto,item_lock);
           } else {
             LOG("Registry", pERROR)
               << "Null TH2F with key = " << name << " - not copied";
           }
     } else if (var_type == "P4TTree") {
           TTree * tree = registry.GetTTree(name);
           if(tree) {
               //TTree * ctree = new TTree(*tree);
               TTree * ctree = tree->CopyTree("1");
               LOG("Registry", pDEBUG) << ctree->GetName();
               cri = new RegistryItem<TTree*>(ctree,item_lock);
           } else {
             LOG("Registry", pERROR)
               << "Null TTree with key = " << name << " - not copied";
           }

     } else {}

     pair<string, RegistryItemI *> reg_entry(name, cri);
     fRegistry.insert(reg_entry);
   }

  fName             = registry.fName;
  fIsReadOnly       = registry.fIsReadOnly;
  fInhibitItemLocks = registry.fInhibitItemLocks;
}
//____________________________________________________________________________
void Registry::Init(void)
{
// initialize registry properties

  fName              = "unnamed";
  fIsReadOnly        = false;
  fInhibitItemLocks  = false;
}
//____________________________________________________________________________
void Registry::Clear(bool force)
{
// clean all registry entries

  if(!force) {
   if(this->IsLocked()) {
      LOG("Registry", pWARN) << "Registry is locked. Can't clear its entries";
      return;
    }
  }

  map<string, RegistryItemI *>::iterator reg_iter;
  for(reg_iter = fRegistry.begin(); reg_iter != fRegistry.end(); reg_iter++) {
     string          name = reg_iter->first;
     RegistryItemI * item = reg_iter->second;

     delete item;
     item = 0;
  }
  fRegistry.clear();
}
//____________________________________________________________________________

