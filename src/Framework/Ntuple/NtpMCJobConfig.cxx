//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TFolder.h>
#include <TObjString.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCJobConfig.h"
#include "Framework/Registry/Registry.h"
#include "Framework/Utils/StringUtils.h"

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

  LOG("Ntp", pNOTICE)
        << "Converting configuration registries to TFolders";

  fConfig = gROOT->GetRootFolder()->AddFolder("gconfig","GENIE configs");
  gROOT->GetListOfBrowsables()->Add(fConfig,"gconfig");

  AlgConfigPool * algconf = AlgConfigPool::Instance();

  const vector<string> & vconfkeys = algconf->ConfigKeyList();
  vector<string>::const_iterator keyiter;

  for(keyiter = vconfkeys.begin(); keyiter != vconfkeys.end(); ++keyiter) {

    string key = *keyiter;

    LOG("Ntp",pDEBUG) << "Current configuration registry key" << key;

    vector<string> vkey = utils::str::Split(key,"/");
    assert(vkey.size()==2);
    string alg_name  = vkey[0];
    string param_set = vkey[1];

    LOG("Ntp",pDEBUG)
         << "alg_name: " << alg_name << ", param_set: " << param_set;

    if( !(fConfig->FindObject(alg_name.c_str())) ) {
      LOG("Ntp",pDEBUG) << "Adding new folder for alg: " << alg_name;
      fConfig->AddFolder(alg_name.c_str(), "");
    }
    TFolder * alg_folder = (TFolder *) fConfig->FindObject(alg_name.c_str());

    LOG("Ntp",pDEBUG) << "Adding folder for param set: " << param_set;
    TFolder * config_folder = alg_folder->AddFolder(param_set.c_str(), "");

    LOG("Ntp",pDEBUG) << "Accessing Registry & converting it to TFolder";
    Registry * config_registry = algconf->FindRegistry(key);
    config_registry->CopyToFolder(config_folder);
  }

  return fConfig;
}
//____________________________________________________________________________
