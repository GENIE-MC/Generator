//____________________________________________________________________________
/*!

\class    genie::RegistryItemI

\brief    Registry item pABC

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

\created  May 04, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_I_H_
#define _REGISTRY_ITEM_I_H_

#include <iostream>

#include "Framework/Registry/RegistryItemTypeId.h"

using std::ostream;

namespace genie {

class RegistryItemI
{
public:
  virtual ~RegistryItemI() { }

  virtual RegistryItemI * Clone    (void)      const = 0;
  virtual RgType_t        TypeInfo (void)      const = 0;
  virtual bool            IsLocked (void)      const = 0;
  virtual void            Lock     (void)            = 0;
  virtual void            UnLock   (void)            = 0;
  virtual bool            IsLocal  (void)      const = 0;
  virtual void            SetLocal (bool)            = 0;
  virtual void            Print    (ostream &) const = 0;

protected:

  RegistryItemI() { }
};

}      // genie namespace

#endif // _REGISTRY_ITEM_I_H_
