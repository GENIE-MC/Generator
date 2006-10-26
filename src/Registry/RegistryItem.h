//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_H_
#define _REGISTRY_ITEM_H_

#include <map>
#include <string>
#include <ostream>

#include "Registry/RegistryItemI.h"
#include "Registry/RegistryItemTypeDef.h"

namespace genie {

template<class T> class RegistryItem;
template<class T>
       ostream & operator << (ostream & stream, const RegistryItem<T> & rec);

template<class T> class RegistryItem : public RegistryItemI {

public:
  RegistryItem() { };
  RegistryItem(T item, bool locked = false);
  RegistryItem(const RegistryItem * ri);
  ~RegistryItem();

  RgType_t  TypeInfo (void) const;
  const T & Data     (void) const { return fItem;         }
  void      Lock     (void)       { fIsLocked = true;     }
  void      UnLock   (void)       { fIsLocked = false;    }
  bool      IsLocked (void) const { return fIsLocked;     }

  void Print(ostream& stream) const;

  friend ostream & operator <<
                        <T>(ostream & stream, const RegistryItem<T> & rec);

private:

  T    fItem;
  bool fIsLocked;
};

}      // genie namespace

#endif // _REGISTRY_ITEM_H_
