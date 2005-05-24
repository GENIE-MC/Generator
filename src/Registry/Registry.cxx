//____________________________________________________________________________
/*!

\class    genie::Registry

\brief    Concrete implementation of the RegistryI interface.
          Provides the container for algorithm configuration parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004
 
*/
//____________________________________________________________________________

#include <cassert>
#include <iomanip>

#include "Messenger/Messenger.h"
#include "Registry.h"

using namespace genie;

using std::pair;
using std::setw;
using std::setfill;
using std::istream;
using std::endl;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const Registry & registry) 
 {
   registry.Print(stream);

   return stream;
 }
}
//____________________________________________________________________________
Registry::Registry() :
fName              ( string("unnamed") ),
fIsReadOnly        ( false             ),
fInhibitItemLocks  ( false             )
{

}
//____________________________________________________________________________
Registry::Registry(const char * name, bool isReadOnly) :
fName             ( string(name) ),
fIsReadOnly       ( isReadOnly   ),
fInhibitItemLocks ( false        )
{

}
//____________________________________________________________________________
Registry::Registry(const Registry & registry)
{
  this->Copy(registry);
}
//____________________________________________________________________________
Registry::~Registry()
{
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
void Registry::operator() (const char * key, int item)
{
  Set(string(key), item);
}
//____________________________________________________________________________
void Registry::operator() (string key, int item)
{
  Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (const char * key, bool item)
{
  Set(string(key), item);
}
//____________________________________________________________________________
void Registry::operator() (string key, bool item)
{
  Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (const char * key, double item)
{
  Set(string(key), item);
}
//____________________________________________________________________________
void Registry::operator() (string key, double item)
{
  Set(key, item);
}
//____________________________________________________________________________
void Registry::operator() (const char * key, string item)
{
  Set(string(key), item);
}
//____________________________________________________________________________
void Registry::operator() (const char * key, const char * item)
{
  Set(string(key), string(item));
}
//____________________________________________________________________________
void Registry::operator() (string key, string item)
{
  Set(key, item);
}
//____________________________________________________________________________
Registry * Registry::operator = (const Registry & reg)
{
  return new Registry(reg);
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
bool Registry::ItemIsLocked(const char * key) const
{
  return this->ItemIsLocked( string(key) );
}
//____________________________________________________________________________
void Registry::LockItem(string key)
{
  if( this->Exists(key) ) {
    
     map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
     entry->second->Lock();

  } else { 
     LOG("Registry", pWARN)
               << "\n ********** Can't lock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
void Registry::LockItem(const char * key)
{
  this->LockItem( string(key) );
}
//____________________________________________________________________________
void Registry::UnLockItem(string key)
{
  if( this->Exists(key) ) {
    
     map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);
     entry->second->UnLock();

  } else {
    LOG("Registry", pWARN)
             << "\n ********** Can't unlock non-existem item [" << key << "]";
  }
}
//____________________________________________________________________________
void Registry::UnLockItem(const char * key)
{
  this->UnLockItem( string(key) );
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
bool Registry::CanSetItem(const char * key) const
{
  return this->CanSetItem( string(key) );
}
//____________________________________________________________________________
void Registry::Set(const char * key, bool item)
{
  this->Set(string(key),item);
}
//____________________________________________________________________________
void Registry::Set(string key, bool item)
{
  if( this->CanSetItem(key) ) {

    bool lock = this->ItemIsLocked(key); // store, could be true but inhibited

    this->DeleteEntry(key);

    RegistryItem<bool> * reg_item = new RegistryItem<bool>(item,lock);

    pair<string, RegistryItemI *> config_entry(key, reg_item);

    fRegistry.insert(config_entry);

  } else {
     LOG("Registry", pWARN)
              << "\n ********** Registry item [" << key << "] can not be set";
  }
}
//____________________________________________________________________________
void Registry::Set(const char * key, int item)
{
  this->Set(string(key),item);
}
//____________________________________________________________________________
void Registry::Set(string key, int item)
{
  if( this->CanSetItem(key) ) {

    bool lock = this->ItemIsLocked(key); // store, could be true but inhibited

    this->DeleteEntry(key);

    RegistryItem<int> * reg_item = new RegistryItem<int>(item,lock);

    pair<string, RegistryItemI *> config_entry(key, reg_item);

    fRegistry.insert(config_entry);

  } else {
     LOG("Registry", pWARN)
              << "\n ********** Registry item [" << key << "] can not be set";
  }
}
//____________________________________________________________________________
void Registry::Set(const char * key, double item)
{
  this->Set(string(key),item);
}
//____________________________________________________________________________
void Registry::Set(string key, double item)
{
  if( this->CanSetItem(key) ) {

    bool lock = this->ItemIsLocked(key); // store, could be true but inhibited

    this->DeleteEntry(key);

    RegistryItem<double> * reg_item = new RegistryItem<double>(item,lock);

    pair<string, RegistryItemI *> config_entry(key, reg_item);

    fRegistry.insert(config_entry);

  } else {
     LOG("Registry", pWARN)
              << "\n ********** Registry item [" << key << "] can not be set";
  }   
}
//____________________________________________________________________________
void Registry::Set(const char * key, const char * item)
{
  this->Set(string(key), string(item));
}
//____________________________________________________________________________
void Registry::Set(const char * key, string item)
{
  this->Set(string(key), item);
}
//____________________________________________________________________________
void Registry::Set(string key, string item)
{
  if( this->CanSetItem(key) ) {
  
    bool lock = this->ItemIsLocked(key); // store, could be true but inhibited

    this->DeleteEntry(key);

    RegistryItem<string> * reg_item = new RegistryItem<string>(item,lock);

    pair<string, RegistryItemI *> config_entry(key, reg_item);

    fRegistry.insert(config_entry);

  } else {
     LOG("Registry", pWARN)
              << "\n ********** Registry item [" << key << "] can not be set";
  }   
}
//____________________________________________________________________________
void Registry::Get(const char * key, bool & item) const
{
  this->Get(string(key),item);
}
//____________________________________________________________________________
void Registry::Get(string key, bool & item) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);

  RegistryItem<bool> * ri = dynamic_cast<RegistryItem<bool>*> (entry->second);

  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(const char * key, int & item) const
{
  this->Get(string(key),item);
}
//____________________________________________________________________________
void Registry::Get(string key, int & item) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);

  RegistryItem<int> * ri = dynamic_cast< RegistryItem<int> * > (entry->second);

  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(const char * key, double & item) const
{
  this->Get(string(key),item);
}
//____________________________________________________________________________
void Registry::Get(string key, double & item) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);

  RegistryItem<double> * ri =
                          dynamic_cast<RegistryItem<double>*> (entry->second);

  item = ri->Data();
}
//____________________________________________________________________________
void Registry::Get(const char * key, string & item) const
{
  this->Get(string(key),item);
}
//____________________________________________________________________________
void Registry::Get(string key, string & item) const
{
  map<string, RegistryItemI *>::const_iterator entry = fRegistry.find(key);

  RegistryItem<string> * ri =
                          dynamic_cast<RegistryItem<string>*> (entry->second);

  item = ri->Data();
}
//____________________________________________________________________________
bool Registry::GetBool(string key) const
{
  bool value;
 
  this->Get(key, value);

  return value;
}
//____________________________________________________________________________
bool Registry::GetBool(const char * key) const
{
  return this->GetBool( string(key) );
}
//____________________________________________________________________________
int Registry::GetInt(string key) const
{
  int value;
 
  this->Get(key, value);

  return value;
}
//____________________________________________________________________________
int Registry::GetInt(const char * key) const
{
  return this->GetInt( string(key) );
}
//____________________________________________________________________________
double Registry::GetDouble(string key) const
{
  double value;
 
  this->Get(key, value);

  return value;
}
//____________________________________________________________________________
double Registry::GetDouble(const char * key) const
{
  return this->GetDouble( string(key) );
}
//____________________________________________________________________________
string Registry::GetString(string key) const
{
  string value;
 
  this->Get(key, value);

  return value;
}
//____________________________________________________________________________
string Registry::GetString(const char * key) const
{
  return this->GetString( string(key) );
}
//____________________________________________________________________________
bool Registry::Exists(const char * key) const
{
  return this->Exists( string(key) );
}
//____________________________________________________________________________
bool Registry::Exists(string key) const
{
  if (fRegistry.count(key) == 1) return true;
  else                           return false;
}
//____________________________________________________________________________
bool Registry::DeleteEntry(const char * key)
{
  return this->DeleteEntry( string(key) );
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
void Registry::SetName(const char * name)
{
  this->SetName( string(name) );
}
//____________________________________________________________________________
void Registry::SetName(string name)
{
  if(! fIsReadOnly) fName = name;
  else {
     LOG("Registry", pWARN)
              << "\n ********** Registry is locked - Can not change its name";
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
           << "\n ********** Key: " << key0
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

   stream << endl;
   
   stream << " |-->[-] number of entries    : "
                       << setfill(' ') << setw(3) << fRegistry.size() << endl;

   map<string, RegistryItemI *>::const_iterator r_iter;
                      
   for(r_iter = fRegistry.begin(); r_iter != fRegistry.end(); r_iter++) {

     string var_name = string("> ") + r_iter->first;
     string var_type = (r_iter->second)->TypeInfo().name();
     
     stream << "      |" << setfill('-') << setw(45) << var_name << " ["
            << setfill(' ') << setw(2) << var_type << "] ";
                                     
     r_iter->second->Print(stream);
   }
}
//____________________________________________________________________________
void Registry::Copy(const Registry & registry)
{
  fName             = registry.fName;
  fIsReadOnly       = registry.fIsReadOnly;
  fInhibitItemLocks = registry.fInhibitItemLocks;

  map<string, RegistryItemI *>::const_iterator reg_iter;

  for(reg_iter = registry.fRegistry.begin();
                          reg_iter != registry.fRegistry.end(); reg_iter++) {

     string name = reg_iter->first;
     RegistryItemI * ri = reg_iter->second;

     string var_type  = ri->TypeInfo().name();
     bool   item_lock = ri->IsLocked();

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
     else {}
     
     pair<string, RegistryItemI *> reg_entry(name, cri);
     
     fRegistry.insert(reg_entry);
   }
}
//____________________________________________________________________________
