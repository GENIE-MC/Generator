//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 12, 2004

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
  string frame(100,'~');
  cout << endl << frame;
  cout << "\nAlgFactory singleton dtor: "
       << "Reporting on algorithms/configurations used during the last job: "
       << endl << frame;
  cout << *this;
  cout << frame << endl;

  cout << "AlgFactory singleton dtor: "
                          << "Deleting all owned algorithmic objects" << endl;
  map<string, Algorithm *>::iterator alg_iter;
  for(alg_iter = fAlgPool.begin(); alg_iter != fAlgPool.end(); ++alg_iter) {
    Algorithm * alg = alg_iter->second;
    if(alg) delete alg;
    alg = 0;
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
const Algorithm * AlgFactory::GetAlgorithm(string alg_name, string param_set)
{
//! Manages the instantiation and "storage/retrieval" of algorithms.
//! These algorithms are owned by the factory and it hands over (to the client)
//! a "const Algorithm *" that can be dynamically casted to the requested
//! Algorithm Interface (eg. XSecAlgorithmI, DecayModelI, PdfModelI, etc...)

  string key = alg_name + "/" + param_set;

  SLOG("AlgFactory", pDEBUG)
                      << "Algorithm: " << key << " requested from AlgFactory";

  if( fAlgPool.count(key) == 1 ) {

     LOG("AlgFactory", pDEBUG) << key << " algorithm found in memory";

     map<string, Algorithm *>::const_iterator alg_iter = fAlgPool.find(key);
     return alg_iter->second;

  } else {
     Algorithm * alg_base = this->InstantiateAlgorithm(alg_name, param_set);

     //-- cache the algorithm for future use
     if(alg_base) {
        pair<string, Algorithm *> key_alg_pair(key, alg_base);
        fAlgPool.insert(key_alg_pair);
     } else {
        LOG("AlgFactory", pFATAL)
                        "Algorithm " << key << " could not be instantiated";
        exit(1);
     }
     return alg_base;
  }
  return 0;
}
//____________________________________________________________________________
Algorithm * AlgFactory::AdoptAlgorithm(string alg_name, string param_set) const
{
//! Hands over an algorithm instance that is owned by the client.
//! The client can alter this object (eg. reconfigure) but the AlgFactory does
//! not keep track of it and the client is responsible for deleting it.

   Algorithm * alg_base = InstantiateAlgorithm(alg_name, param_set);
   return alg_base;
}
//____________________________________________________________________________
Algorithm * AlgFactory::InstantiateAlgorithm(
                                      string alg_name, string param_set) const
{
//! Instantiate the requested object based on the registration of its TClass
//! through the generated ROOT dictionaries
//! The class of any object instantiated here must have a LinkDef entry.

  LOG("AlgFactory", pDEBUG) << "Instantiating algorithm = " << alg_name;

  // Get object through gROOT->GetClass() and cast it to the Algorithm base
  // class (ABC)
  TClass * tclass = gROOT->GetClass(alg_name.c_str());
  if(!tclass) {
     LOG("AlgFactory", pERROR) 
             << "Failed instantiating algorithm = " << alg_name;
     return 0;
  }
  void * vd_base = tclass->New();
  Algorithm * alg_base = (Algorithm *) (vd_base);

  LOG("AlgFactory", pDEBUG) << "Setting Configuration Set = " << param_set;

  // Set the configuration Registry that corresponds to the input param_set
  if(strcmp(param_set.c_str(),"NoConfig")!=0) alg_base->Configure(param_set);

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

