//____________________________________________________________________________
/*!

\class    genie::AlgFactory

\brief    Algorithm Factory.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2004
*/
//____________________________________________________________________________

#ifndef _ALG_FACTORY_H_
#define _ALG_FACTORY_H_

#include <map>
#include <string>

using std::map;
using std::pair;
using std::string;

namespace genie {

class Algorithm;

class AlgFactory {

public:

  static AlgFactory * Instance();

  const Algorithm * GetAlgorithm(
                         string alg_name, string param_set="NoConfig");
                         
  Algorithm * AdoptAlgorithm(
                          string alg_name, string param_set="NoConfig") const;
  
private:

  AlgFactory(); 
  AlgFactory(const AlgFactory & alg_factory);
  virtual ~AlgFactory();

  Algorithm * InstantiateAlgorithm(string alg_name, string param_set) const;

  static AlgFactory * fInstance;

  map<string, Algorithm *> fAlgPool; ///< alg_name/param_set -> Algorithm 
  
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (AlgFactory::fInstance !=0) {
            delete AlgFactory::fInstance;
            AlgFactory::fInstance = 0;
         }
      }
  };

  friend struct Cleaner;
};

}      // genie namespace

#endif // _ALG_FACTORY_H_
