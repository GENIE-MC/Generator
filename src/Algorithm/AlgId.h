//____________________________________________________________________________
/*!

\class    genie::AlgId

\brief    Algorithm ID (algorithm name + configuration set name)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 20, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _ALGID_H_
#define _ALGID_H_

#include <string>
#include <iostream>

#include "Registry/RegistryTypesDef.h"

using std::string;
using std::ostream;

namespace genie {

class AlgId {

public:

  AlgId();
  AlgId(string name, string config);
  AlgId(const AlgId & id);
  AlgId(const RgAlg & registry_item);
 ~AlgId();

  string Name   (void) const { return fName;   }
  string Config (void) const { return fConfig; }
  string Key    (void) const { return fKey;    }

  void   SetId     (string name, string config="");
  void   SetName   (string name);
  void   SetConfig (string config);
  void   Copy      (const AlgId & id);
  void   Copy      (const RgAlg & registry_item);
  void   Print     (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const AlgId & alg);

private:

  void Init      (void);
  void UpdateKey (void);

  string fName;   ///< Algorithm name (including namespaces)
  string fConfig; ///< Configuration set name
  string fKey;    ///< Unique key: namespace::alg_name/alg_config
};

}       // genie namespace

#endif  // _ALGID_H_
