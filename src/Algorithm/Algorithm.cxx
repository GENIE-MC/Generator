//____________________________________________________________________________
/*!

\class    genie::Algorithm

\brief    Algorithm abstract base class. 

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004
 
*/
//____________________________________________________________________________

#include "Algorithm/AlgFactory.h"
#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const Algorithm & alg)
  {
     // print algorithm name & parameter-set
     stream << alg.fName << "/" << alg.fParamSet << endl;

     // print algorithm configuration
     stream << *(alg.fConfig);

     return stream;
  }
}
//____________________________________________________________________________
Algorithm::Algorithm()
{
  fConfigIsOwned = false; 
  fConfig        = 0;
} 
//____________________________________________________________________________
Algorithm::Algorithm(const char * param_set) :
fParamSet(param_set)
{ 
  fConfigIsOwned = false; 
  fConfig        = 0;
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
void Algorithm::Configure(string param_set)
{
// Configure the Algorithm looking up at the ConfigPool singleton for a
// configuration Registry corresponding to the input named parameter set.

  fParamSet = param_set;

  FindConfig();
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
                      << this->Name() << "/" << this->ParamSet() 
                                                    << " at the ConfigPool";
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
const Algorithm * Algorithm::SubAlg(string alg_key, string config_key) const
{
// Returns the sub-algorithm pointed to this algorithm's XML config file using
// the the values of the alg/config keys
// This method asserts the existence of these keys in the XML config.

  LOG("Algorithm", pDEBUG)
               << "Fetching sub-algorithm within algorithm: " << this->Name();
  LOG("Algorithm", pDEBUG)
              << "Asserting existence of keys: ["
                                    << alg_key << "], [" << config_key << "]";

  fConfig->AssertExistence(alg_key, config_key);

  string alg_name  = fConfig->GetString(alg_key);
  string param_set = fConfig->GetString(config_key);

  LOG("Algorithm", pDEBUG)
              << "Keys retrieved these values: alg="
                                     << alg_name << ", config=" << param_set;
                                     
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
               << "Fetching sub-algorithm within algorithm: " << this->Name();
  LOG("Algorithm", pDEBUG)
              << "Testing existence of keys: ["
                                    << alg_key << "], [" << config_key << "]";

  bool keys_exist = (fConfig->Exists(alg_key) && fConfig->Exists(config_key));

  if(keys_exist) {
     string alg_name  = fConfig->GetString(alg_key);
     string param_set = fConfig->GetString(config_key);

     LOG("Algorithm", pDEBUG)
              << "Keys retrieved these values: alg="
                                     << alg_name << ", config=" << param_set;
                                     
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

