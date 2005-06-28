//____________________________________________________________________________
/*!

\class    genie::AlgConfigPool

\brief    A singleton class which holds all configuration registries assembled
          from XML configuration files for all (agorithm-name, parameter-set)
          pairs. Any algorithmic object can get an instance of the algorithm
          config. pool and query it to learn its configuration parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#include <TSystem.h>

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "XML/Xml2ConfigFileList.h"
#include "XML/Xml2Registry.h"

using namespace genie;

using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<<(ostream & stream, const AlgConfigPool & config_pool)
  {
    config_pool.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
AlgConfigPool * AlgConfigPool::fInstance = 0;
//____________________________________________________________________________
AlgConfigPool::AlgConfigPool()
{
  if( ! LoadXMLConfig() )

  LOG("AlgConfigPool", pERROR) << "Could not load XML config file";

  fInstance =  0;
}
//____________________________________________________________________________
AlgConfigPool::~AlgConfigPool()
{
  fInstance = 0;
}
//____________________________________________________________________________
AlgConfigPool * AlgConfigPool::Instance()
{
  if(fInstance == 0) {

    static AlgConfigPool::Cleaner cleaner;

    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new AlgConfigPool;
  }

  return fInstance;
}
//____________________________________________________________________________
bool AlgConfigPool::LoadXMLConfig(void)
{
  LOG("AlgConfigPool", pINFO)
        << "AlgConfigPool late initialization: Loading all XML config. files";

  //-- get base GENIE directory from $GENIE environmental variable

  string base_dir = string( gSystem->Getenv("GENIE") );

  //-- build base config directory and master config filename

  string config_dir    = base_dir + string("/config");
  string master_config = base_dir + string("/config/master_conf.xml");

  //-- check XML structure of MASTER_CONFIG file

  if(Xml2ConfigFileList::CheckXmlParsing(master_config)
                                                   != kXmlOK) return false;

  //-- build a "algorithm" -> "XML config file" look-up map

  map<string, string> conf_files =
                           Xml2ConfigFileList::ReadFileList(master_config);

  //-- loop over all XML config files and read all named configuration
  //   sets for each algorithm

  map<string, string>::const_iterator conf_file_iter;

  for(conf_file_iter = conf_files.begin();
                    conf_file_iter != conf_files.end(); ++conf_file_iter) {

    string alg_name    = conf_file_iter->first;
    string config_file = config_dir + "/" + conf_file_iter->second;

    LOG("AlgConfigPool", pINFO) << alg_name << " ---> " << config_file;

    map<string, Registry *> alg_configs =
                                 Xml2Registry::ReadAlgConfigs(config_file);

    map<string, Registry *>::const_iterator conf_iter;
    for(conf_iter = alg_configs.begin();
                             conf_iter != alg_configs.end(); ++conf_iter) {

        string     alg_param_set = alg_name + "/" + conf_iter->first;
        Registry * conf_registry = conf_iter->second;

        conf_registry->SetName(alg_param_set);
        conf_registry->Lock();

        pair<string, Registry *> alg_conf(alg_param_set, conf_registry);
        fRegistryPool.insert(alg_conf);
    }
  }

  return true;
};
//____________________________________________________________________________
Registry * AlgConfigPool::FindRegistry(const Algorithm * algorithm) const
{
  return FindRegistry( algorithm->Name(), algorithm->ParamSet() );
}
//____________________________________________________________________________
Registry* AlgConfigPool::FindRegistry(string alg_name, string param_set) const
{
  string key = alg_name + "/" + param_set;

  LOG("AlgConfigPool", pDEBUG) << "Searching for registry with key " << key;

  if( fRegistryPool.count(key) == 1 ) {

     map<string, Registry *>::const_iterator config_entry =
                                                   fRegistryPool.find(key);
     return config_entry->second;

  } else {
     LOG("AlgConfigPool", pWARN) << "No config registry for key " << key;
     return 0;
  }
}
//____________________________________________________________________________
void AlgConfigPool::Print(ostream & stream) const
{
  stream << "\n[-] ALGORITHM CONFIGURATION POOL: ";

  typedef map<string, Registry *>::const_iterator  sregIter;
  typedef map<string, Registry *>::size_type       sregSize;

  sregSize size = fRegistryPool.size();

  stream << size << " configuration sets found" << endl;
  for(sregIter iter = fRegistryPool.begin();
                                       iter != fRegistryPool.end(); iter++) {
     stream << iter->first     << endl;
     stream << *(iter->second) << endl;
  }
}
//____________________________________________________________________________

