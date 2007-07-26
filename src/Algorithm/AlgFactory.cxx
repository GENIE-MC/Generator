//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <iostream>
#include <cstdlib>

#include <TROOT.h>
#include <TClass.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/Algorithm.h"
#include "Messenger/Messenger.h"

using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
  ostream & operator<<(ostream & stream, const AlgFactory & algf)
  {
    algf.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
AlgFactory * AlgFactory::fInstance = 0;
//____________________________________________________________________________
AlgFactory::AlgFactory()
{
  fInstance =  0;
}
//____________________________________________________________________________
AlgFactory::~AlgFactory()
{
  // clean-up 
  cout << "AlgFactory singleton dtor: "
       << "Deleting all owned algorithmic objects used at this job" << endl;

  map<string, Algorithm *>::iterator alg_iter;
  for(alg_iter = fAlgPool.begin(); alg_iter != fAlgPool.end(); ++alg_iter) {
    Algorithm * alg = alg_iter->second;
    if(alg) {
      cout << "- Deleting algorithm: " << alg->Id() << endl;
      delete alg;
      alg = 0;
    }
  }
  fAlgPool.clear();
  fInstance = 0;
}
//____________________________________________________________________________
AlgFactory * AlgFactory::Instance()
{
  if(fInstance == 0) {
    static AlgFactory::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new AlgFactory;
  }
  return fInstance;
}
//____________________________________________________________________________
const Algorithm * AlgFactory::GetAlgorithm(const AlgId & algid)
{
//! Manages the instantiation and "storage/retrieval" of algorithms.
//! These algorithms are owned by the factory and it hands over (to the client)
//! a "const Algorithm *" that can be dynamically casted to the requested
//! Algorithm Interface (eg. XSecAlgorithmI, DecayModelI, PdfModelI, etc...)

  return this->GetAlgorithm(algid.Name(), algid.Config());
}
//____________________________________________________________________________
const Algorithm * AlgFactory::GetAlgorithm(string name, string config)
{
  string key = name + "/" + config;

  SLOG("AlgFactory", pDEBUG)
      << "Algorithm: " << key << " requested from AlgFactory";

  map<string, Algorithm *>::const_iterator alg_iter = fAlgPool.find(key);
  bool found = (alg_iter != fAlgPool.end());

  if(found) {
     LOG("AlgFactory", pDEBUG) << key << " algorithm found in memory";
     return alg_iter->second;
  } else {
     //-- instantiate the factory
     Algorithm * alg_base = this->InstantiateAlgorithm(name,config);

     //-- cache the algorithm for future use
     if(alg_base) {
        pair<string, Algorithm *> key_alg_pair(key, alg_base);
        fAlgPool.insert(key_alg_pair);
     } else {
        LOG("AlgFactory", pFATAL)
            << "Algorithm: " << key << " could not be instantiated";
        exit(1);
     }
     return alg_base;
  }
  return 0;
}
//____________________________________________________________________________
Algorithm * AlgFactory::AdoptAlgorithm(const AlgId & algid) const
{
//! Hands over an algorithm instance that is owned by the client.
//! The client can alter this object (eg. reconfigure) but the AlgFactory does
//! not keep track of it and the client is responsible for deleting it.

  return this->AdoptAlgorithm(algid.Name(), algid.Config());
}
//____________________________________________________________________________
Algorithm * AlgFactory::AdoptAlgorithm(string name, string config) const
{
   Algorithm * alg_base = InstantiateAlgorithm(name, config);
   return alg_base;
}
//____________________________________________________________________________
void AlgFactory::ForceReconfiguration(void)
{
  LOG("AlgFactory", pINFO) 
       << "Forcing reconfiguration of all owned algorithms";

  map<string, Algorithm *>::iterator alg_iter = fAlgPool.begin();
  for( ; alg_iter != fAlgPool.end(); ++alg_iter) {
    Algorithm * alg = alg_iter->second;
    string config = alg->Id().Config();
    bool skip_conf = (config=="NoConfig" || config=="");
    if(!skip_conf) alg->Configure(config);
  }
}
//____________________________________________________________________________
Algorithm * AlgFactory::InstantiateAlgorithm(string name, string config) const
{
//! Instantiate the requested object based on the registration of its TClass
//! through the generated ROOT dictionaries
//! The class of any object instantiated here must have a LinkDef entry.

  // Get object through ROOT's TROOT::GetClass() mechanism
  LOG("AlgFactory", pDEBUG) << "Instantiating algorithm = " << name;

  TClass * tclass = gROOT->GetClass(name.c_str());
  if(!tclass) {
     LOG("AlgFactory", pERROR) 
         << "Failed instantiating algorithm = " << name;
     return 0;
  }
  void * vd_base = tclass->New();
  Algorithm * alg_base = (Algorithm *) (vd_base);

  // Configure the instantiated algorithm

  LOG("AlgFactory", pDEBUG) << "Setting Configuration Set = " << config;

  bool skip_conf = (config=="NoConfig" || config=="");
  if(skip_conf) {
    LOG("AlgFactory", pDEBUG) << "Skipping algorithm configuration step!";
  } 
  else alg_base->Configure(config);

  return alg_base;
}
//____________________________________________________________________________
void AlgFactory::Print(ostream & stream) const
{
  string frame(100,'.');

  stream << endl;
  map<string, Algorithm *>::const_iterator alg_iter;
  for(alg_iter = fAlgPool.begin(); alg_iter != fAlgPool.end(); ++alg_iter) {
    const Algorithm * alg = alg_iter->second;
    stream << frame << endl;
    stream << "Used algorithm: " << alg->Id() << endl;
    stream << "Printing config:";
    stream << alg->GetConfig();
  }

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry & gc = *(confp->GlobalParameterList());

  stream << frame << endl;
  stream << "Printing global parameters list:";
  stream << gc;
}
//____________________________________________________________________________

