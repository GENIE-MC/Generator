//____________________________________________________________________________
/*!

\class    genie::AlgFactory

\brief    Algorithm Factory.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2004
*/
//____________________________________________________________________________

#include <TROOT.h>
#include <TClass.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/Algorithm.h"
#include "Messenger/Messenger.h"

namespace genie {

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
  
     Algorithm * alg_base = InstantiateAlgorithm(alg_name, param_set);

     //-- cache the algorithm for future use
     if(alg_base) {
        pair<string, Algorithm *> key_alg_pair(key, alg_base);
        
        fAlgPool.insert(key_alg_pair);
        
     } else {

        LOG("AlgFactory", pFATAL)
                               << key << " could not be instantiated" << ENDL;

        assert(false);
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

  void * base = gROOT->GetClass(alg_name.c_str())->New();

  // Cast to the Algorithm Abstract Base Class (ABC)

  Algorithm * alg_base = (Algorithm *) base;

  // Set the configuration Registry that corresonds to the input param_set

  LOG("AlgFactory", pDEBUG) << "Setting Configuration Set = " << param_set;

  if( strcmp(param_set.c_str(),"NoConfig") != 0)
                                     alg_base->Configure(param_set);

  return alg_base;
}
//____________________________________________________________________________

} // genie namespace 
