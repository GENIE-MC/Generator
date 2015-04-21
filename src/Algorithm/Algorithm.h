//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 02, 2004

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_

#include <string>
#include <iostream>
#include <cassert>
#include <map>

#include "Algorithm/AlgStatus.h"
#include "Algorithm/AlgCmp.h"
#include "Algorithm/AlgId.h"
#include "Registry/Registry.h"
#include "Registry/RegistryItemTypeDef.h"

using std::string;
using std::ostream;
using std::map;

namespace genie {

class Algorithm;
typedef map <string, Algorithm *>                 AlgMap;
typedef map <string, Algorithm *>::iterator       AlgMapIter;
typedef map <string, Algorithm *>::const_iterator AlgMapConstIter;
typedef pair<string, Algorithm *>                 AlgMapPair;

class Algorithm {

public:
  virtual ~Algorithm();

  //! Configure the algorithm
  virtual void Configure (const Registry & config);  

  //! Configure the algorithm 
  virtual void Configure (string config);            

  //! Lookup configuration from the config pool
  virtual void FindConfig (void);   

  //! Get configuration registry
  virtual const Registry & GetConfig(void) const { return *fConfig; }

  //! Get a writeable version of an owned configuration Registry.
  Registry * GetOwnedConfig(void);

  //! Get algorithm ID
  virtual const AlgId & Id(void) const { return fID; }

  //! Get algorithm status
  virtual AlgStatus_t GetStatus(void) const { return fStatus; }

  //! Allow reconfigration after initializaton?
  //! Algorithms may opt-out, if reconfiguration is not necessary, 
  //! to improve event reweighting speed.
  virtual bool AllowReconfig(void) const { return fAllowReconfig; }

  //! Compare with input algorithm
  virtual AlgCmp_t Compare(const Algorithm * alg) const;

  //! Set algorithm ID
  virtual void SetId(const AlgId & id);
  virtual void SetId(string name,  string config);

  //! Access the sub-algorithm pointed to by the input key, either from the 
  //! local pool or from AlgFactory's pool
  const Algorithm * SubAlg(const RgKey & registry_key) const;

  //! Clone the configuration registry looked up from the configuration pool
  //! and take its ownership
  void AdoptConfig (void);

  //! Take ownership of the algorithms subtructure (sub-algorithms,...)
  //! by copying them from the AlgFactory pool to the local pool
  //! Also bring all the configuration variables to the top level config Registry.
  //! This can be used to group together a series of algorithms & their
  //! configurations and extract (a clone of) this group from the shared 
  //! pools. Having a series of algorithms/configurations behaving as a
  //! monolithic block, with a single point of configuration (the top level)
  //! is to be used when bits & pieces of GENIE are used in isolation for
  //! data fitting or reweighting
  void AdoptSubstructure (void);

  //! Print algorithm info
  virtual void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const Algorithm & alg);

protected:
  Algorithm();
  Algorithm(string name);
  Algorithm(string name, string config);

  void Initialize         (void);
  void DeleteConfig       (void);
  void DeleteSubstructure (void);

  bool         fAllowReconfig; ///<
  bool         fOwnsConfig;    ///< true if it owns its config. registry
  bool         fOwnsSubstruc;  ///< true if it owns its substructure (sub-algs,...)
  AlgId        fID;            ///< algorithm name and configuration set
  Registry *   fConfig;        ///< config. (either owned or pointing to config pool)
  AlgStatus_t  fStatus;        ///< algorithm execution status
  AlgMap *     fOwnedSubAlgMp; ///< local pool for owned sub-algs (taken out of the factory pool)
};

}       // genie namespace
#endif  // _ALGORITHM_H_
