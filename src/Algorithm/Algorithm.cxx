//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 20, 2009 - CA
   Added fAllowReconfig private data member and AllowReconfig() method.
   Algorithms can set this method to opt-out of reconfiguration. Speeds up 
   reweighting if algorithms (that don't need to be reconfigured) opt out.
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
  this->DeleteConfig();
  this->DeleteSubstructure();
}
//____________________________________________________________________________
void Algorithm::Configure(const Registry & config)
{
// Configure the Algorithm using the input configuration Registry

  LOG("Algorithm", pNOTICE) << "Input configuration: " << config;

  if(!fOwnsConfig) {
    fConfig     = new Registry(config);
    fOwnsConfig = true;
  } else {
    fConfig->Copy(config);
  }
  LOG("Algorithm", pNOTICE) << "Copied configuration: " << *fConfig;

  if(!fOwnsSubstruc) return;              // doesn't own substructure
  if(fOwnedSubAlgMp->size()==0) return;   // no sub-algorithms

  LOG("Algorithm", pNOTICE) << "Configuring algorithms stored at local pool";

  // loop over local pool algorithms
  AlgMapIter alg_iter = fOwnedSubAlgMp->begin();
  for( ; alg_iter != fOwnedSubAlgMp->end(); ++alg_iter) {
    string      alg_key = alg_iter->first;
    Algorithm * alg     = alg_iter->second;
    if(!alg) {
      LOG("Algorithm", pERROR) 
       << "Key: " << alg_key << " points to a null algorithm at local pool";
      continue; 
    }
    LOG("Algorithm", pNOTICE) << "Configuring alg: " << alg->Id().Key();

    // get local pool algorithm's owned config
    Registry r(alg->GetConfig());

    const RgIMap & rgmap = fConfig->GetItemMap();
    RgIMapConstIter reg_iter = rgmap.begin();
    for( ; reg_iter != rgmap.end(); ++reg_iter) {
       RgKey reg_key = reg_iter->first;
       if(reg_key.find(alg_key+"/") == string::npos) continue;

       int i0 = reg_key.find_first_of("/")+1;
       int i1 = reg_key.length();
       RgKey new_reg_key = reg_key.substr(i0,i1);
       LOG("Algorithm", pDEBUG) 
 	 << "Item with key: " << reg_key 
                   << " is copied 1-level down with key: " << new_reg_key;

       RegistryItemI * ri = reg_iter->second;
       RgIMapPair key_item_pair(new_reg_key, ri->Clone());
       r.Set(key_item_pair);
    }
    alg->Configure(r);
  }
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
  Registry * config = pool->FindRegistry(this);

  if(!config)
     // notify & keep whatever config Registry was used before.
      LOG("Algorithm", pWARN)
                   << "No Configuration available for "
                               << this->Id().Key() << " at the ConfigPool";
  else {
     // check if its already owns a configuration Registry & delete it;
     if(fOwnsConfig && fConfig) delete fConfig;

     // point to the configuration Registry retrieved from the ConfigPool
     // and raise the "not-owned" flag.
     fConfig        = config;
     fOwnsConfig = false;

     LOG("Algorithm", pDEBUG) << "\n" << *fConfig;
  }
}
//____________________________________________________________________________
Registry * Algorithm::GetOwnedConfig(void)
{
  if(!fOwnsConfig) return 0;
  return fConfig;
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
  stream << "\nAlgorithm Key: "  << this->fID.Key();
  stream << " - Owns Config: "   << ((fOwnsConfig)   ? "[true]" : "[false]");
  stream << " - Owns Substruc: " << ((fOwnsSubstruc) ? "[true]" : "[false]");

  // print algorithm configuration
  const Registry & r = this->GetConfig();
  stream << r;

  if(fOwnsSubstruc) {
    AlgMapConstIter iter = fOwnedSubAlgMp->begin();
    for(; iter!=fOwnedSubAlgMp->end(); ++iter) {
      Algorithm * alg = iter->second;
      stream << "<Next algorithm is owned by : " << this->fID.Key() << ">";
      stream << *alg;
    }
  }
}
//____________________________________________________________________________
void Algorithm::Initialize(void)
{
// Algorithm initialization
//
  fAllowReconfig  = true;
  fOwnsConfig     = false;
  fOwnsSubstruc   = false;
  fConfig         = 0;
  fOwnedSubAlgMp  = 0;
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
//
  LOG("Algorithm", pINFO)
    << "Fetching sub-alg within alg: " << this->Id().Key()
                                   << " pointed to by key: " << registry_key;

  //-- if the algorithm owns its substructure:
  //      return the sub-algorithm from the local pool
  //
  if(fOwnsSubstruc) {
     AlgMapConstIter iter = fOwnedSubAlgMp->find(registry_key);
     if(iter!=fOwnedSubAlgMp->end()) return iter->second;
     LOG("Algorithm", pERROR) 
       << "Owned sub-alg pointed to by key: " << registry_key 
                       << " was not found within alg: " << this->Id().Key();
     return 0;
  }

  //-- if the algorithm dowes not own its substructure:
  //      return the sub-algorithm from the AlgFactory's pool
  //
  fConfig->AssertExistence(registry_key);

  // retrieve the algorithm item corresponding to key
  RgAlg alg = fConfig->GetAlg(registry_key);
  LOG("Algorithm", pINFO)
     << "Registry key: " << registry_key << " points to algorithm: " << alg;

  // retrieve the Algorithm object from the the Algorithm factory
  AlgFactory * algf = AlgFactory::Instance();
  const Algorithm * algbase = algf->GetAlgorithm(alg.name, alg.config);
  assert(algbase);

  return algbase;
}
//____________________________________________________________________________
void Algorithm::AdoptConfig(void)
{
  LOG("Algorithm", pNOTICE)
    << this->Id().Key() << " is taking ownership of its configuration";

  if(fOwnsConfig) {
    LOG("Algorithm", pWARN)
      << this->Id().Key() << " already owns its configuration!";
    return;
  }
  Registry & configuration = *fConfig;
  this->Configure(configuration);
}
//____________________________________________________________________________
void Algorithm::AdoptSubstructure(void)
{
// Take ownership of the algorithms subtructure (sub-algorithms,..) by copying 
// them from the AlgFactory pool to the local pool. Also bring all the 
// configuration variables to the top level. See the header for more details.
// A secial naming convention is required for configuration parameter keys 
// for parameters belonging to sub-algorithms (or sub-algorithms of these 
// sub-algorithms and so on...). 
// The convention is: "sub-alg-key/sub-sub-alg-key/.../original name"
// This is a recursive method: The AdoptSubtructure()of all sub-algorithms is
// invoked.
//
  LOG("Algorithm", pNOTICE) 
     << "Algorithm: " << this->Id().Key() << " is adopting its substructure";

  Registry deep_config;
  deep_config.UnLock();
  deep_config.SetName(this->Id().Key());

  //  deep_config.SetName(this->Id().Config() + ";D");
  //  fID.SetConfig(this->Id().Config() + ";D");

  if(fOwnsSubstruc) this->DeleteSubstructure();

  fOwnedSubAlgMp = new AlgMap;
  fOwnsSubstruc  = true;

  AlgFactory * algf = AlgFactory::Instance();

  const RgIMap & rgmap = fConfig->GetItemMap();

  RgIMapConstIter iter = rgmap.begin();
  for( ; iter != rgmap.end(); ++iter) {

    RgKey reg_key = iter->first;
    RegistryItemI * ri = iter->second;

    if(ri->TypeInfo() == kRgAlg) {
        LOG("Algorithm", pDEBUG) 
                 << "Found sub-algorithm pointed to by " << reg_key;
        RgAlg reg_alg = fConfig->GetAlg(reg_key);
        AlgId id(reg_alg);

        LOG("Algorithm", pDEBUG) << "Adopting sub-algorithm = " << id.Key();
        Algorithm * subalg = algf->AdoptAlgorithm(id.Name(),id.Config());
        subalg->AdoptSubstructure();

        LOG("Algorithm", pDEBUG) << "Adding sub-algorithm to local pool";
        AlgMapPair key_alg_pair(reg_key, subalg);
        fOwnedSubAlgMp->insert(key_alg_pair);

        LOG("Algorithm", pDEBUG) << "Appending its config";
        const Registry & r = subalg->GetConfig();
        RgKey prefix = reg_key + "/";
        deep_config.Append(r,prefix);

    } else {
       LOG("Algorithm", pDEBUG) << "Adding parameter with key = " << reg_key;
       RgIMapPair key_item_pair(reg_key, ri->Clone());
       deep_config.Set(key_item_pair);
    }
  }
  this->Configure(deep_config);
}
//____________________________________________________________________________
void Algorithm::DeleteConfig(void)
{
  // there is nothing to delete if the configuration is not owned but is 
  // rather looked up from the configuration pool
  //
  if(!fOwnsConfig) return;

  // delete owned configuration registry
  //
  if(fConfig) {
    delete fConfig; 
    fConfig=0;
  }
}
//____________________________________________________________________________
void Algorithm::DeleteSubstructure(void)
{
  // there is nothing to delete if the sub-algorithms are not owned but rather
  // taken from the AlgFactory's pool
  //
  if(!fOwnsSubstruc) return;

  // delete local algorithm pool
  //
  AlgMapIter iter = fOwnedSubAlgMp->begin();
  for( ; iter != fOwnedSubAlgMp->end(); ++iter) {
    Algorithm * alg = iter->second;
    if(alg) {
      delete alg;
      alg=0;
    }
  }
  delete fOwnedSubAlgMp;
  fOwnedSubAlgMp = 0;
}
//____________________________________________________________________________
