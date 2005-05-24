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

#ifndef _REGISTRY_H_
#define _REGISTRY_H_

#include <map>
#include <string>
#include <iostream>

#include "RegistryItem.h"

using std::map;
using std::string;
using std::ostream;

namespace genie {

class Registry {

public:

  Registry();
  Registry(const char * name, bool isReadOnly = true);
  Registry(const Registry &);
  virtual ~Registry();

  friend ostream & operator << (ostream & stream, const Registry & registry);
  
  Registry * operator = (const Registry & reg);

  void operator () (const char * key,  int          item);
  void operator () (string       key,  int          item);
  void operator () (const char * key,  bool         item);
  void operator () (string       key,  bool         item);
  void operator () (const char * key,  double       item);  
  void operator () (string       key,  double       item);
  void operator () (const char * key,  string       item);  
  void operator () (string       key,  string       item);
  void operator () (const char * key,  const char * item);
  
  //! Registry locks

  void   Lock               (void);                   ///< locks the registry
  void   UnLock             (void);                   ///< unlocks the registry (doesn't unlock items)
  bool   IsLocked           (void) const;             ///< checks registry lock
  
  //! Registry item locks

  void   InhibitItemLocks   (void);                   ///< override individual item locks
  void   RestoreItemLocks   (void);                   ///< restore individual item locks
  bool   ItemLocksAreActive (void) const;             ///< check if individual item locks are active
  void   LockItem           (string       key);       ///< locks the registry item
  void   LockItem           (const char * key);       ///< locks the registry item
  void   UnLockItem         (string       key);       ///< unlocks the registry item
  void   UnLockItem         (const char * key);       ///< unlocks the registry item
  bool   ItemIsLocked       (string       key) const; ///< check item lock
  bool   ItemIsLocked       (const char * key) const; ///< check item lock

  //! Methods to set/retrieve Registry values
  
  void   Set         (const char * key, bool         item);
  void   Set         (string       key, bool         item);
  void   Set         (const char * key, int          item);
  void   Set         (string       key, int          item);
  void   Set         (const char * key, double       item);
  void   Set         (string       key, double       item);
  void   Set         (const char * key, const char * item);
  void   Set         (const char * key, string       item);
  void   Set         (string       key, string       item);
  void   Get         (const char * key, bool &       item) const;
  void   Get         (string       key, bool &       item) const;
  void   Get         (const char * key, int &        item) const;
  void   Get         (string       key, int &        item) const;
  void   Get         (const char * key, double &     item) const;
  void   Get         (string       key, double &     item) const;
  void   Get         (const char * key, string &     item) const;
  void   Get         (string       key, string &     item) const;
  bool   GetBool     (string       key) const;
  bool   GetBool     (const char * key) const;
  int    GetInt      (string       key) const;
  int    GetInt      (const char * key) const;
  double GetDouble   (string       key) const;
  double GetDouble   (const char * key) const;
  string GetString   (string       key) const;
  string GetString   (const char * key) const;

  int    NEntries    (void) const;
  bool   Exists      (const char * key) const;
  bool   Exists      (string       key) const;
  bool   CanSetItem  (string       key) const;
  bool   CanSetItem  (const char * key) const;
  bool   DeleteEntry (const char * key);
  bool   DeleteEntry (string       key);
  void   SetName     (const char * name);
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

}        // namespace

#endif   // _REGISTRY_H_
