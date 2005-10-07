//____________________________________________________________________________
/*!

\class   genie::NtpMCJobConfig

\brief   Stores the GENIE configuration in ROOT TFolders along with the
         output event tree

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <cassert>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TFolder.h>
#include <TObjString.h>

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCJobConfig.h"
#include "Registry/Registry.h"
#include "Utils/StringUtils.h"

using std::vector;
using std::string;
using namespace genie;

ClassImp(NtpMCJobConfig)

//____________________________________________________________________________
NtpMCJobConfig::NtpMCJobConfig()
{
  fConfig = 0;
}
//____________________________________________________________________________
NtpMCJobConfig::~NtpMCJobConfig()
{

}
//____________________________________________________________________________
TFolder * NtpMCJobConfig::Load(void)
{
  if (fConfig) delete fConfig;
  fConfig = 0;

  LOG("NtpMCConf", pNOTICE)
                << "Converting configuration registries to TFolders";

  fConfig = gROOT->GetRootFolder()->AddFolder("gconfig","GENIE configs");
  gROOT->GetListOfBrowsables()->Add(fConfig,"gconfig");

  AlgConfigPool * algconf = AlgConfigPool::Instance();

  const vector<string> & vconfkeys = algconf->ConfigKeyList();
  vector<string>::const_iterator keyiter;

  for(keyiter = vconfkeys.begin(); keyiter != vconfkeys.end(); ++keyiter) {

    string key = *keyiter;

    LOG("NtpMCConf",pDEBUG) << "Current configuration registry key" << key;

    vector<string> vkey = utils::str::Split(key,"/");
    assert(vkey.size()==2);
    string alg_name  = vkey[0];
    string param_set = vkey[1];

    LOG("NtpMCConf",pDEBUG)
                << "alg_name: " << alg_name << ", param_set: " << param_set;

    if( !(fConfig->FindObject(alg_name.c_str())) ) {
      LOG("NtpMCConf",pDEBUG) << "Adding new folder for alg: " << alg_name;
      fConfig->AddFolder(alg_name.c_str(), "");
    }
    TFolder * alg_folder = (TFolder *) fConfig->FindObject(alg_name.c_str());

    LOG("NtpMCConf",pDEBUG) << "Adding folder for param set: " << param_set;
    TFolder * config_folder = alg_folder->AddFolder(param_set.c_str(), "");

    LOG("NtpMCConf",pDEBUG) << "Accessing Registry & converting it to TFolder";
    Registry * config_registry = algconf->FindRegistry(key);
    config_registry->CopyToFolder(config_folder);
  }

  return fConfig;
}
//____________________________________________________________________________
