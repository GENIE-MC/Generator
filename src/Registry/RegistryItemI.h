//____________________________________________________________________________
/*!

\class    genie::RegistryItemI

\brief    Registry item pABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_I_H_
#define _REGISTRY_ITEM_I_H_

//#include <typeinfo>
#include <iostream>

#include "Registry/RegistryItemTypeId.h"

using std::ostream;
//using std::type_info;

namespace genie {

class RegistryItemI 
{
public:
  virtual ~RegistryItemI() { }

  //virtual const type_info & TypeInfo (void) const = 0;

  virtual RegistryItemI * Clone    (void)      const = 0;
  virtual RgType_t        TypeInfo (void)      const = 0;
  virtual void            Lock     (void)            = 0;
  virtual void            UnLock   (void)            = 0;
  virtual bool            IsLocked (void)      const = 0;
  virtual void            Print    (ostream &) const = 0;

protected:

  RegistryItemI() { }
};

}      // genie namespace

#endif // _REGISTRY_ITEM_I_H_
