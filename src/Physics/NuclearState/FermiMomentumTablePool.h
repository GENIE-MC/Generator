//____________________________________________________________________________
/*!

\class    genie::FermiMomentumTablePool

\brief    Singleton class to load & serve tables of Fermi momentum constants

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 18, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FERMI_MOMENTUM_TABLE_POOL_H_
#define _FERMI_MOMENTUM_TABLE_POOL_H_

#include <map>
#include <string>

#include "Framework/Conventions/XmlParserStatus.h"

using std::map;
using std::string;

namespace genie {

class FermiMomentumTable;

class FermiMomentumTablePool
{
public:
  static FermiMomentumTablePool * Instance (void);

  const FermiMomentumTable * GetTable (string name);

private:
  FermiMomentumTablePool();
  FermiMomentumTablePool(const FermiMomentumTablePool & fmc);
  virtual ~FermiMomentumTablePool();

  bool LoadTables (void);
  XmlParserStatus_t ParseXMLTables (string filename);

  static FermiMomentumTablePool * fInstance;

  map<string, FermiMomentumTable *> fKFSets; // kf tables

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (FermiMomentumTablePool::fInstance !=0) {
            delete FermiMomentumTablePool::fInstance;
            FermiMomentumTablePool::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _FERMI_MOMENTUM_TABLE_POOL_H_
