//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldUserData

\brief    A singleton that manages the DBTable<T> objects, extracted from the
          NuVld RDBMS using its GUI, and makes them available to other objects
          such as the NuVldKernel.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 21, 2005
*/
//_____________________________________________________________________________

#ifndef _NUVLD_USER_DATA_H_
#define _NUVLD_USER_DATA_H_

#include <vector>
#include <string>

#include <TObject.h>

#include "DBUtils/DBTable.h"
#include "DBUtils/DBTableType.h"
#include "DBUtils/DBTableStack.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "Facades/NGCardPairList.h"

using std::vector;
using std::string;

namespace genie {
namespace nuvld {

class NuVldUserData : public TObject {

public:

   static NuVldUserData * Instance();

   NGCardPair *                      NeuGenCardPair      (void) { return fCurrNeuGenCardPair;  }
   NGCardPairList *                  NeuGenCardPairStack (void) { return fNeuGenCardPairStack; }
   DBTableStack<vXSecTableRow> *     NuXSecStack         (void) { return fNuXSecTableStack;    }
   DBTableStack<eDiffXSecTableRow> * ElDiffXSecStack     (void) { return fElDiffXSecTableStack;}

   bool CurrDBTableIsNull (void) const { return (fCurrDBTable == 0); }

   void SetCurrDBTable (genie::nuvld::DBTable<vXSecTableRow> *     table);
   void SetCurrDBTable (genie::nuvld::DBTable<eDiffXSecTableRow> * table);
      
   void DelCurrDBTable (void);
   
   void AddCurrDBTableToStack   (string name);
   void SetStackedDBTableAsCurr (string name);
   void DelStackedDBTable       (string name);

   DBTableType_t CurrDBTableType(void) const { return fCurrDBTableType; }
   DBTableBase * CurrDBTable    (void) const { return fCurrDBTable;     }
   
   DBTable<vXSecTableRow> *     NuXSec     (void);
   DBTable<eDiffXSecTableRow> * ElDiffXSec (void);

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
   NGCardPair *  fCurrNeuGenCardPair;
   
   NGCardPairList *                  fNeuGenCardPairStack;
   DBTableStack<vXSecTableRow> *     fNuXSecTableStack;
   DBTableStack<eDiffXSecTableRow> * fElDiffXSecTableStack;
      
ClassDef(NuVldUserData, 0)
};

} // nuvld namespace
} // genie namespace

#endif

