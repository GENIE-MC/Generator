//____________________________________________________________________________
/*!

\class    genie::Registry

\brief    A registry. Provides the container for algorithm configuration
          parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _REGISTRY_H_
#define _REGISTRY_H_

#include <map>
#include <string>
#include <iostream>

#include "Registry/RegistryItem.h"

class TH1F;
class TH2F;
class TTree;

using std::map;
using std::pair;
using std::string;
using std::ostream;

namespace genie {

class Registry {

public:

  Registry();
  Registry(string name, bool isReadOnly = true);
  Registry(const Registry &);
  virtual ~Registry();

  friend ostream & operator << (ostream & stream, const Registry & registry);

  Registry * operator = (const Registry & reg);

  void operator () (string       key,  int          item);
  void operator () (string       key,  bool         item);
  void operator () (string       key,  double       item);
  void operator () (string       key,  const char * item);
  void operator () (string       key,  string       item);

  //! Registry locks

  void   Lock      (void);       ///< locks the registry
  void   UnLock    (void);       ///< unlocks the registry (doesn't unlock items)
  bool   IsLocked  (void) const; ///< checks registry lock

  //! Registry item locks

  void   InhibitItemLocks   (void);             ///< override individual item locks
  void   RestoreItemLocks   (void);             ///< restore individual item locks
  bool   ItemLocksAreActive (void) const;       ///< check if item locks are active
  void   LockItem           (string key);       ///< locks the registry item
  void   UnLockItem         (string key);       ///< unlocks the registry item
  bool   ItemIsLocked       (string key) const; ///< check item lock

  //! Methods to set/retrieve Registry values

  void   Set (pair<string, genie::RegistryItemI *> entry);
  void   Get (string key, const RegistryItemI * item) const;

  void   Set (string key, bool        item);
  void   Set (string key, int         item);
  void   Set (string key, double      item);
  void   Set (string key, string      item);
  void   Set (string key, const char* item);
  void   Set (string key, TH1F *      item);
  void   Set (string key, TH2F *      item);
  void   Set (string key, TTree *     item);
  void   Get (string key, bool &      item) const;
  void   Get (string key, int &       item) const;
  void   Get (string key, double &    item) const;
  void   Get (string key, string &    item) const;
  void   Get (string key, TH1F *      item) const;
  void   Get (string key, TH2F *      item) const;
  void   Get (string key, TTree *     item) const;

  bool   GetBool     (string key) const;
  int    GetInt      (string key) const;
  double GetDouble   (string key) const;
  string GetString   (string key) const;
  TH1F * GetTH1F     (string key) const;
  TH2F * GetTH2F     (string key) const;
  TTree* GetTTree    (string key) const;

  int    NEntries    (void) const;
  bool   Exists      (string key) const;
  bool   CanSetItem  (string key) const;
  bool   DeleteEntry (string key);
  void   SetName     (string name);
  string Name        (void) const;
  void   Print       (ostream & stream) const;
  void   Copy        (const Registry &);

  //! Assert the existence or registry items

  void   AssertExistence (string key0) const;
  void   AssertExistence (string key0, string key1) const;
  void   AssertExistence (string key0, string key1, string key2) const;

private:

  string fName;
  bool   fIsReadOnly;
  bool   fInhibitItemLocks;

  map<string, RegistryItemI *> fRegistry;
};

template<class T> void SetRegistryItem(Registry * r, string key, T   item);

}        // namespace

#endif   // _REGISTRY_H_
