//____________________________________________________________________________
/*!

\class    genie::AlgConfigPool

\brief    A singleton class holding all configuration registries built while
          parsing all loaded XML configuration files. 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALG_CONFIG_POOL_H_
#define _ALG_CONFIG_POOL_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Registry/Registry.h"

using std::map;
using std::vector;
using std::string;
using std::ostream;

namespace genie {

class AlgConfigPool;
ostream & operator << (ostream & stream, const AlgConfigPool & cp);

class AlgConfigPool {

public:
  static AlgConfigPool * Instance();

  Registry * FindRegistry (string key)                        const;
  Registry * FindRegistry (string alg_name, string param_set) const;
  Registry * FindRegistry (const Algorithm * algorithm)       const;
  Registry * FindRegistry (const AlgId & algid)               const;

  Registry * GlobalParameterList(void) const;
  Registry * CommonList( const string & file_id, const string & set_name ) const;
  Registry * TuneGeneratorList(void) const;
 

  const vector<string> & ConfigKeyList (void) const;

  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const AlgConfigPool & cp);

private:
  AlgConfigPool();
  AlgConfigPool(const AlgConfigPool & config_pool);
  virtual ~AlgConfigPool();

  // methods for loading all algorithm XML configuration files
  string BuildConfigKey      (string alg_name, string param_set) const;
  string BuildConfigKey      (const Algorithm * algorithm) const;
  bool   LoadAlgConfig       (void);
  bool   LoadMasterConfig    (void);
  bool   LoadGlobalParamLists(void);
  bool   LoadCommonLists( const string & file_id );
  bool   LoadTuneGeneratorList(void);
  bool   LoadSingleAlgConfig (string alg_name, string file_name);
  bool   LoadRegistries      (string key_base, string file_name, string root);
  void   AddConfigParameter  (Registry * r, string pt, string pn, string pv);
  void   AddBasicParameter   (Registry * r, string pt, string pn, string pv);
  void   AddRootObjParameter (Registry * r, string pt, string pn, string pv);

  static AlgConfigPool * fInstance;

  map<string, Registry *> fRegistryPool;  ///< algorithm/param_set -> Registry
  map<string, string>     fConfigFiles;   ///< algorithm -> XML config file
  vector<string>          fConfigKeyList; ///< list of all available configuration keys
  string                  fMasterConfig;  ///< lists config files for all algorithms

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (AlgConfigPool::fInstance !=0) {
            delete AlgConfigPool::fInstance;
            AlgConfigPool::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _ALG_CONFIG_POOL_H_
