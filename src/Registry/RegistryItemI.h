//____________________________________________________________________________
/*!

\class    genie::RegistryItemI

\brief    Registry Item pABC (interface)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#ifndef _REGISTRY_ITEM_I_H_
#define _REGISTRY_ITEM_I_H_

#include <typeinfo>
#include <iostream>

using std::ostream;
using std::type_info;

namespace genie {

class RegistryItemI 
{
public:

  virtual ~RegistryItemI() { }

  virtual const type_info & TypeInfo (void)             const = 0;
  virtual void              Lock     (void)                   = 0;
  virtual void              UnLock   (void)                   = 0;
  virtual bool              IsLocked (void)             const = 0;
  virtual void              Print    (ostream & stream) const = 0;

protected:

  RegistryItemI()          { }
};

}      // genie namespace

#endif // _REGISTRY_ITEM_I_H_
