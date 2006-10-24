//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 02, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <vector>
#include <string>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

using std::vector;
using std::string;
using std::endl;

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const Algorithm & alg)
  {
    alg.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
Algorithm::Algorithm()
{
  this->Initialize();
}
//____________________________________________________________________________
Algorithm::Algorithm(string name)
{
  this->Initialize();
  fID.SetId(name);
}
//____________________________________________________________________________
Algorithm::Algorithm(string name, string config)
{
  this->Initialize();
  fID.SetId(name,config);
  this->FindConfig();
}
//____________________________________________________________________________
Algorithm::~Algorithm()
{
  //if(fConfigIsOwned && fConfig) delete fConfig; -- check
}
//____________________________________________________________________________
void Algorithm::Configure(const Registry & config)
{
// Configure the Algorithm using the input configuration Registry

  // check if its already owns a configuration Registry & delete it

  if(fConfigIsOwned && fConfig) {

    LOG("Algorithm", pINFO) << "Deleting old owned configuration";

    //delete fConfig; -- check
  }

  // create a new configuration Registry & raise the "owned" flag

  fConfig        = new Registry(config);
  fConfigIsOwned = true;

  LOG("Algorithm", pINFO) << *fConfig;
}
//____________________________________________________________________________
void Algorithm::Configure(string config)
{
// Configure the Algorithm looking up at the ConfigPool singleton for a
// configuration Registry corresponding to the input named parameter set.

  fID.SetConfig(config);

  this->FindConfig();
}
//____________________________________________________________________________
void Algorithm::FindConfig(void)
{
// Finds its configration Registry from the ConfigPool and gets a pointer to
// it. If the Registry comes from the ConfigPool then the Algorithm does not
// own its configuration (the ConfigPool singleton keeps the ownership and the
// responsibility to -eventually- delete all the Registries it instantiates
// by parsing the XML config files).

  AlgConfigPool * pool = AlgConfigPool::Instance();

  Registry * config = pool->FindRegistry( this );

  if(!config)
        // notify & keep whatever config Registry was used before.
        LOG("Algorithm", pWARN)
                   << "No Configuration available for "
                               << this->Id().Key() << " at the ConfigPool";
  else {
        // check if its already owns a configuration Registry & delete it;
        if(fConfigIsOwned && fConfig) delete fConfig;

        // point to the configuration Registry retrieved from the ConfigPool
        // and raise the "not-owned" flag.
        fConfig        = config;
        fConfigIsOwned = false;

        LOG("Algorithm", pDEBUG) << ENDL << *fConfig;
  }
}
//____________________________________________________________________________
AlgCmp_t Algorithm::Compare(const Algorithm * algo) const
{
// Compares itself with the input algorithm

  string alg1    = this->Id().Name();
  string config1 = this->Id().Config();
  string alg2    = algo->Id().Name();
  string config2 = algo->Id().Config();

  if(alg1 == alg2)
  {
    if(config1 == config2) return kAlgCmpIdentical;
    else                   return kAlgCmpDiffConfig;
  }
  else return kAlgCmpDiffAlg;

  return kAlgCmpUnknown;
}
//____________________________________________________________________________
void Algorithm::SetId(const AlgId & id)
{
  fID.Copy(id);
}
//____________________________________________________________________________
void Algorithm::SetId(string name, string config)
{
  fID.SetId(name, config);
}
//____________________________________________________________________________
void Algorithm::Print(ostream & stream) const
{
  // print algorithm name & parameter-set
  stream << this->fID.Key() << endl;

  // print algorithm configuration
  stream << *(this->fConfig);
}
//____________________________________________________________________________
void Algorithm::Initialize()
{
  fConfigIsOwned = false;
  fConfig        = 0;
}
//____________________________________________________________________________
const Algorithm * Algorithm::SubAlg(const RgKey & registry_key) const
{
// Returns the sub-algorithm pointed to this algorithm's XML config file using
// the the values of the key.
// This method asserts the existence of these keys in the XML config.
// Note: Since only 1 parameter is used, the key value should contain both the
// algorithm name and its configuration set according to the usual scheme:
// namespace::algorithm_name/configuration_set

  LOG("Algorithm", pDEBUG)
          << "Fetching sub-algorithm within algorithm: " << this->Id().Key();
  LOG("Algorithm", pDEBUG)
               << "Asserting existence of key: [" << registry_key << "]";

  // assert key existence 
  fConfig->AssertExistence(registry_key);

  // retrieve the algorithm item corresponding to key
  RgAlg alg = fConfig->GetAlg(registry_key);
  LOG("Algorithm", pDEBUG)
     << "Registry key: " << registry_key << " points to algorithm: " << alg;

  //vector<string> algv = str::Split(alg,"/");
  //assert(algv.size()==2);
  //string alg_name  = algv[0];
  //string param_set = algv[1];
  //LOG("Algorithm", pDEBUG) << "Key value split, alg-name  : " << alg_name;
  //LOG("Algorithm", pDEBUG) << "Key value split, param_set : " << param_set;

  // retrieve the Algorithm object from the the Algorithm factory
  AlgFactory * algf = AlgFactory::Instance();
  const Algorithm * algbase = algf->GetAlgorithm(alg.name, alg.config);
  assert(algbase);

  return algbase;
}
//____________________________________________________________________________
/*
const Algorithm * Algorithm::SubAlg(string alg_key, string config_key) const
{
// Returns the sub-algorithm pointed to this algorithm's XML config file using
// the the values of the alg/config keys
// This method asserts the existence of these keys in the XML config.

  LOG("Algorithm", pDEBUG)
          << "Fetching sub-algorithm within algorithm: " << this->Id().Key();
  LOG("Algorithm", pDEBUG)
              << "Asserting existence of keys: ["
                                   << alg_key << "], [" << config_key << "]";

  fConfig->AssertExistence(alg_key, config_key);

  string alg_name  = fConfig->GetString(alg_key);
  string param_set = fConfig->GetString(config_key);

  LOG("Algorithm", pDEBUG)
          << "Input key [" << alg_key << "] points to value: " << alg_name;
  LOG("Algorithm", pDEBUG)
      << "Input key [" << config_key << "] points to value: " << param_set;

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);
  assert(algbase);

  return algbase;
}
//____________________________________________________________________________
const Algorithm * Algorithm::SubAlgWithDefault(string alg_key,
         string config_key, string def_alg_name, string def_config_name) const
{
// Returns the sub-algorithm pointed to this algorithm's XML config file using
// the the values of the alg/config keys
// This method does not assert the existence of these keys in the XML config.
// If the keys are not found it will attempt to load the specified 'Default'
// preconfigured sub-algorithm

  AlgFactory * algf = AlgFactory::Instance();

  LOG("Algorithm", pDEBUG)
          << "Fetching sub-algorithm within algorithm: " << this->Id().Key();
  LOG("Algorithm", pDEBUG)
              << "Testing existence of keys: ["
                                   << alg_key << "], [" << config_key << "]";

  bool keys_exist = (fConfig->Exists(alg_key) && fConfig->Exists(config_key));

  if(keys_exist) {
     string alg_name  = fConfig->GetString(alg_key);
     string param_set = fConfig->GetString(config_key);

     LOG("Algorithm", pDEBUG)
              << "Input keys point to algorithm: "
                                        << alg_name << "/" << param_set;

     const Algorithm* algbase = algf->GetAlgorithm(alg_name, param_set);
     assert(algbase);

     return algbase;
  }

  LOG("Algorithm", pDEBUG)
       << "Keys were not found. Loading defaults";
  LOG("Algorithm", pDEBUG)
       << "Defaults: alg=" << def_alg_name << ", config=" << def_config_name;

  const Algorithm * defalgbase =
                           algf->GetAlgorithm(def_alg_name, def_config_name);
  assert(defalgbase);

  return defalgbase;
}
//____________________________________________________________________________
*/
