//____________________________________________________________________________
/*!

\class    genie::AlgFactory

\brief    The GENIE Algorithm Factory.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 12, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _ALG_FACTORY_H_
#define _ALG_FACTORY_H_

#include <map>
#include <string>
#include <iostream>

#include "Framework/Algorithm/AlgId.h"

using std::map;
using std::pair;
using std::string;
using std::ostream;

namespace genie {

class Algorithm;
class AlgFactory;

ostream & operator << (ostream & stream, const genie::AlgFactory & algf);

class AlgFactory {

public:
  static AlgFactory * Instance();

  //! Instantiates, configures and returns a pointer to the specified algorithm.
  //! The algorithm is placed at the AlgFactory pool (is owned by the factory)
  //! from where it will be looked up at subsequent calls.
  const Algorithm * GetAlgorithm (const AlgId & algid);
  const Algorithm * GetAlgorithm (string name, string conf="Default");

  //! Like GetAlgorithm() but the algorithm is not placed at the AlgFactory pool
  //! and its ownership is transfered to the caller
  Algorithm * AdoptAlgorithm (const AlgId & algid) const;
  Algorithm * AdoptAlgorithm (string name, string conf="Default") const;

  //! Forces a reconfiguration of all algorithms kept at the factory pool.
  //! The algorithms look up their nominal configuration from the config pool.
  //! Use that to propagate modifications made directly at the config pool.
  void ForceReconfiguration(bool ignore_alg_opt_out=false);

  //! print algorithm factory
  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const AlgFactory & algf);

private:
  AlgFactory();
  AlgFactory(const AlgFactory & alg_factory);
  virtual ~AlgFactory();

  //! method instantiating (based on TClass * TROOT::GetClass(name)) 
  //! & configuring algorithmic objects
  Algorithm * InstantiateAlgorithm(string name, string config) const;

  //! sinleton's self
  static AlgFactory * fInstance;

  //! 'algorithm key' (namespace::name/config) -> 'algorithmic object' map
  map<string, Algorithm *> fAlgPool;

  //! singleton cleaner
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
