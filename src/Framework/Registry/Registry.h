//____________________________________________________________________________
/*!

\class    genie::Registry

\brief    A registry. Provides the container for algorithm configuration
          parameters.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REGISTRY_H_
#define _REGISTRY_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "Framework/Registry/RegistryItem.h"
#include "Framework/Registry/RegistryItemTypeDef.h"

class TH1F;
class TH2F;
class TTree;
class TFolder;

using std::map;
using std::vector;
using std::pair;
using std::string;
using std::ostream;

namespace genie {

// Type definitions
//
typedef map <RgKey, RegistryItemI *>                 RgIMap;
typedef pair<RgKey, RegistryItemI *>                 RgIMapPair;
typedef map <RgKey, RegistryItemI *>::size_type      RgIMapSizeType;
typedef map <RgKey, RegistryItemI *>::iterator       RgIMapIter;
typedef map <RgKey, RegistryItemI *>::const_iterator RgIMapConstIter;
typedef vector<RgKey>                                RgKeyList;

// Templated utility methods to set/get registry items
//
class Registry;
template<class T> 
  void SetRegistryItem(Registry * r, RgKey key, T item);
template<class T> 
  T GetValueOrUseDefault(
     Registry * r, RgKey key, T def, bool set_def=true);

//
//
ostream & operator << (ostream & stream, const Registry & registry);

class Registry {

public:
  // Ctor's & dtor
  //
  Registry();
  Registry(string name, bool isReadOnly = true);
  Registry(const Registry &);
  virtual ~Registry();

  // Overloaded registry operators (<<, (), = , +=)
  //
  friend ostream & operator << (ostream & stream, const Registry & registry);

  Registry & operator =  (const Registry & reg);
  Registry & operator += (const Registry & reg);

  void operator () (RgKey key,  int          item);
  void operator () (RgKey key,  bool         item);
  void operator () (RgKey key,  double       item);
  void operator () (RgKey key,  const char * item);
  void operator () (RgKey key,  string       item);

  // Registry & registry item locks
  //
  void   Lock                (void);            ///< locks the registry
  void   UnLock              (void);            ///< unlocks the registry (doesn't unlock items)
  bool   IsLocked            (void) const;      ///< checks registry lock
  void   InhibitItemLocks    (void);            ///< override individual item locks
  void   RestoreItemLocks    (void);            ///< restore individual item locks
  bool   ItemLocksAreActive  (void) const;      ///< check if item locks are active
  void   LockItem            (RgKey key);       ///< locks the registry item
  void   UnLockItem          (RgKey key);       ///< unlocks the registry item
  bool   ItemIsLocked        (RgKey key) const; ///< check item lock
  bool   ItemIsLocal         (RgKey key) const; ///< local or global?
  void   OverrideGlobalDef   (RgKey key);       ///< let item override global default   (i.e. a 'local'  item)
  void   LinkToGlobalDef     (RgKey key);       ///< link its value to a global default (i.e. a 'global' item)

  // Methods to set/retrieve Registry values
  //
  void   Set (RgIMapPair entry);
  void   Set (RgKey key, RgBool  item);
  void   Set (RgKey key, RgInt   item);
  void   Set (RgKey key, RgDbl   item);
  void   Set (RgKey key, RgStr   item);
  void   Set (RgKey key, RgAlg   item);
  void   Set (RgKey key, RgCChAr item);
  void   Set (RgKey key, RgH1F   item);
  void   Set (RgKey key, RgH2F   item);
  void   Set (RgKey key, RgTree  item);

  void   Get (RgKey key, const RegistryItemI * & item) const;
  void   Get (RgKey key, RgBool & item) const;
  void   Get (RgKey key, RgInt &  item) const;
  void   Get (RgKey key, RgDbl &  item) const;
  void   Get (RgKey key, RgStr &  item) const;
  void   Get (RgKey key, RgAlg &  item) const;
  void   Get (RgKey key, RgH1F &  item) const;
  void   Get (RgKey key, RgH2F &  item) const;
  void   Get (RgKey key, RgTree & item) const;

  RgBool GetBool      (RgKey key) const;
  RgInt  GetInt       (RgKey key) const;
  RgDbl  GetDouble    (RgKey key) const;
  RgStr  GetString    (RgKey key) const;
  RgAlg  GetAlg       (RgKey key) const;
  RgH1F  GetH1F       (RgKey key) const;
  RgH2F  GetH2F       (RgKey key) const;
  RgTree GetTree      (RgKey key) const;

  RgBool GetBoolDef   (RgKey key, RgBool def_opt, bool set_def=true);
  RgInt  GetIntDef    (RgKey key, RgInt  def_opt, bool set_def=true);
  RgDbl  GetDoubleDef (RgKey key, RgDbl  def_opt, bool set_def=true);
  RgStr  GetStringDef (RgKey key, RgStr  def_opt, bool set_def=true);
  RgAlg  GetAlgDef    (RgKey key, RgAlg  def_opt, bool set_def=true);
  
  RgIMapConstIter SafeFind  (RgKey key) const;

  int    NEntries     (void) const;                     ///< get number of items
  bool   Exists       (RgKey key) const;                ///< item with input key exists?
  bool   CanSetItem   (RgKey key) const;                ///< can I set the specifed item?
  bool   DeleteEntry  (RgKey key);                      ///< delete the spcified item
  void   SetName      (string name);                    ///< set the registry name
  string Name         (void) const;                     ///< get the registry name
  void   Print        (ostream & stream) const;         ///< print the registry to stream
  void   Copy         (const Registry &);               ///< copy the input registry
  void   Append       (const Registry &, RgKey pfx=""); ///< append the input registry. Entries already in the registry are not updated
  void   Merge        (const Registry &, RgKey pfx=""); ///< append the input registry. Entries already in the registry are updated
  void   Clear        (bool force = false);             ///< clear the registry
  void   Init         (void);                           ///< initialize the registry

  RgType_t  ItemType (RgKey key)      const;  ///< return item type
  RgKeyList FindKeys (RgKey key_part) const;  ///< create list with all keys containing 'key_part'

  // Access key->item map
  //
  const RgIMap & GetItemMap(void) const { return fRegistry; }

  // Convert to TFolder (this is the primary mechanism for saving the
  // GENIE configuration in a ROOT file, along with its generated events)
  //
  void CopyToFolder (TFolder * folder) const;

  // Assert the existence or registry items
  //
  void AssertExistence (RgKey key0) const;
  void AssertExistence (RgKey key0, RgKey key1) const;
  void AssertExistence (RgKey key0, RgKey key1, RgKey key2) const;

private:

  RegistryItemI * CloneRegistryItem( const RgKey & key ) const ;   ///< Properly clone a registry Item according to its type

  // Registry's private data members
  //
  string fName;              ///< registry's name
  bool   fIsReadOnly;        ///< is read only?
  bool   fInhibitItemLocks;  ///<
  RgIMap fRegistry;          ///< 'key' -> 'value' map
};

}        // genie namespace

#endif   // _REGISTRY_H_
