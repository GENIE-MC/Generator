//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_

#include <string>
#include <iostream>
#include <cassert>

#include "Algorithm/AlgStatus.h"
#include "Algorithm/AlgCmp.h"
#include "Algorithm/AlgId.h"
#include "Registry/Registry.h"
#include "Registry/RegistryItemTypeDef.h"

using std::string;
using std::ostream;

namespace genie {

class Algorithm {

public:
  virtual ~Algorithm();

  //-- Define the Algorithm interface

  //! Configure the algorithm
  virtual void Configure(const Registry & config);

  //! Configure the algorithm
  virtual void Configure(string config);

  //! Lookup configuration from the algorithm configuration pool
  virtual void FindConfig(void);

  //! Get configuration registry
  virtual const Registry & GetConfig(void) const { return *fConfig; }

  //! Get algorithm ID
  virtual const AlgId & Id(void) const { return fID; }

  //! Get algorithm status
  virtual AlgStatus_t GetStatus(void) const { return fStatus; }

  //! Compare with input algorithm
  virtual AlgCmp_t Compare(const Algorithm * alg) const;

  //! Set algorithm ID
  virtual void SetId(const AlgId & id);

  //! Set algorithm ID
  virtual void SetId(string name,  string config);

  //! Print algorithm info
  virtual void Print(ostream & stream) const;

  //! Print algorithm info
  friend ostream & operator << (ostream & stream, const Algorithm & alg);

protected:
  Algorithm();
  Algorithm(string name);
  Algorithm(string name, string config);

  void Initialize();

  const Algorithm * SubAlg(const RgKey & registry_key) const;
/*
  const Algorithm * SubAlg(const AlgId & aid_id)       const;
  const Algorithm * SubAlgWithDefault(string alg_key, string config_key,
                                string def_alg_name, string def_config_name) const;
*/
  bool         fConfigIsOwned; ///< true if the algorithm owns its config. registry
  AlgId        fID;            ///< algorithm name and configuration set
  Registry *   fConfig;        ///< config. (either owned or pointing to config pool)
  AlgStatus_t  fStatus;        ///< algorithm execution status
};

}       // genie namespace
#endif  // _ALGORITHM_H_
