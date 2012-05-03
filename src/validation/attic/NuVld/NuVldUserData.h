//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldUserData

\brief    A singleton that manages the DBTable<T> objects, extracted from the
          NuVld RDBMS using its GUI, and makes them available to other objects
          such as the NuVldKernel.

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 21, 2005
*/
//_____________________________________________________________________________

#ifndef _NUVLD_USER_DATA_H_
#define _NUVLD_USER_DATA_H_

#include <vector>
#include <string>

#include <TObject.h>

#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBTableType.h"
#include "ValidationTools/NuVld/DBTableStack.h"
#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBSFTableRow.h"

using std::vector;
using std::string;

namespace genie {
namespace nuvld {

class NuVldUserData : public TObject {

public:

   static NuVldUserData * Instance();

   DBTableStack<DBNuXSecTableRow> *     NuXSecStack         (void) { return fNuXSecTableStack;    }
   DBTableStack<DBElDiffXSecTableRow> * ElDiffXSecStack     (void) { return fElDiffXSecTableStack;}
   DBTableStack<DBSFTableRow> *         SFStack             (void) { return fSFTableStack;        }

   bool CurrDBTableIsNull (void) const { return (fCurrDBTable == 0); }

   void SetCurrDBTable (genie::nuvld::DBTable<DBNuXSecTableRow> *     table);
   void SetCurrDBTable (genie::nuvld::DBTable<DBElDiffXSecTableRow> * table);
   void SetCurrDBTable (genie::nuvld::DBTable<DBSFTableRow> *         table);
      
   void DelCurrDBTable (void);
   
   void AddCurrDBTableToStack   (string name);
   void SetStackedDBTableAsCurr (string name);
   void DelStackedDBTable       (string name);

   DBTableType_t CurrDBTableType(void) const { return fCurrDBTableType; }
   DBTableBase * CurrDBTable    (void) const { return fCurrDBTable;     }
   
   DBTable<DBNuXSecTableRow> *     NuXSec     (void);
   DBTable<DBElDiffXSecTableRow> * ElDiffXSec (void);
   DBTable<DBSFTableRow> *         SF         (void);

   unsigned int NStackedDBTables (void) const;
   unsigned int NStackedConfigs  (void) const;

   void ClearStackedConfig(void);
   void ClearStackedDBTables(void);

   const vector<string> * GetStackedDBTableNames(void) const;

   void SaveStack(string root_file) const;
                  
private:

   NuVldUserData();
   NuVldUserData(const NuVldUserData & user_data);
   ~NuVldUserData();

   void Init(void);

   static NuVldUserData * fSelf;

   DBTableType_t fCurrDBTableType;   
   DBTableBase * fCurrDBTable;   
   
   DBTableStack<DBNuXSecTableRow> *     fNuXSecTableStack;
   DBTableStack<DBElDiffXSecTableRow> * fElDiffXSecTableStack;
   DBTableStack<DBSFTableRow> *         fSFTableStack;
      
ClassDef(NuVldUserData, 0)
};

} // nuvld namespace
} // genie namespace

#endif

