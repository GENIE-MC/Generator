//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_H_
#define _REGISTRY_ITEM_H_

#include <map>
#include <string>
#include <ostream>

#include "Framework/Registry/RegistryItemI.h"
#include "Framework/Registry/RegistryItemTypeDef.h"

namespace genie {

template<typename T> class RegistryItem;
template<typename T>
  ostream & operator << (ostream & stream, const RegistryItem<T> & rec);

template<typename T> class RegistryItem : public RegistryItemI {

public:
  RegistryItem() { };
  RegistryItem(T item, bool locked=false, bool local=true);
  RegistryItem(const RegistryItem * ri);
  ~RegistryItem();

  RegistryItemI * Clone    (void) const;
  RgType_t        TypeInfo (void) const;
  const T &       Data     (void) const {  return fItem;       }
  bool            IsLocked (void) const {  return fIsLocked;   }
  void            Lock     (void)       {  fIsLocked = true;   }
  void            UnLock   (void)       {  fIsLocked = false;  }
  bool            IsLocal  (void) const {  return fIsLocal;    }
  void            SetLocal (bool isloc) {  fIsLocal  = isloc;  }

  void Print(ostream& stream) const;

  friend ostream & operator <<
              <T>(ostream & stream, const RegistryItem<T> & rec);

private:

  T    fItem;
  bool fIsLocked;
  bool fIsLocal;
};

}      // genie namespace

#endif // _REGISTRY_ITEM_H_
