//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004

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

using std::string;
using std::ostream;

namespace genie {

class Algorithm {

public:

  virtual ~Algorithm();

  //-- define the Algorithm interface

  virtual void             Configure  (const Registry & config);
  virtual void             Configure  (string config);
  virtual void             FindConfig (void);
  virtual const Registry & GetConfig  (void) const { return *fConfig; }
  virtual const AlgId &    Id         (void) const { return  fID;     }
  virtual AlgStatus_t      GetStatus  (void) const { return  fStatus; }
  virtual AlgCmp_t         Compare    (const Algorithm * alg) const;
  virtual void             SetId      (const AlgId & id);
  virtual void             SetId      (string name,  string config);
  virtual void             Print      (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const Algorithm & alg);

protected:

  Algorithm();
  Algorithm(string name);
  Algorithm(string name, string config);

  void Initialize();

  const Algorithm * SubAlg(string key) const;
  const Algorithm * SubAlg(string alg_key, string config_key) const;
  const Algorithm * SubAlgWithDefault(string alg_key, string config_key,
                                string def_alg_name, string def_config_name) const;

  bool         fConfigIsOwned; ///< true if the algorithm owns its config. registry
  AlgId        fID;            ///< algorithm name and configuration set
  Registry *   fConfig;        ///< config. (either owned or pointing to config pool)
  AlgStatus_t  fStatus;        ///< algorithm execution status
};

}       // genie namespace

#endif  // _ALGORITHM_H_
