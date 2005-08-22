//____________________________________________________________________________
/*!

\class    genie::CharmFractionTablePool

\brief    Singleton class to load & serve a TDatabasePDG.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#ifndef _CHARM_FRACTION_TABLE_POOL_H_
#define _CHARM_FRACTION_TABLE_POOL_H_

#include <map>
#include <string>

#include "Conventions/XmlParserStatus.h"
#include "Fragmentation/CharmFractionTable.h"

using std::map;
using std::string;

namespace genie {

class CharmFractionTablePool
{
public:

  static CharmFractionTablePool * Instance(void);

private:

  CharmFractionTablePool();
  CharmFractionTablePool(const CharmFractionTablePool & config_pool);
  virtual ~CharmFractionTablePool();

  bool              LoadTables     (void);
  XmlParserStatus_t ParseXMLTables (string filename);

  static CharmFractionTablePool * fInstance;

  map<string, CharmFractionTable *> fTablePool;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (CharmFractionTablePool::fInstance !=0) {
            delete CharmFractionTablePool::fInstance;
            CharmFractionTablePool::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _CHARM_FRACTION_TABLE_POOL_H_
