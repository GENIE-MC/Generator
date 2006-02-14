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

#include <TFile.h>

#include "Messenger/Messenger.h"
#include "NuVldGUI/NuVldUserData.h"

using namespace genie::nuvld;

ClassImp(NuVldUserData)

//_____________________________________________________________________________
NuVldUserData * NuVldUserData::fSelf = 0;
//_____________________________________________________________________________
NuVldUserData * NuVldUserData::Instance()
{
  if(fSelf == 0) fSelf = new NuVldUserData;

  return fSelf;     
}
//_____________________________________________________________________________
NuVldUserData::NuVldUserData()
{
  fSelf = 0;

  Init();
}
//_____________________________________________________________________________
void NuVldUserData::Init(void)
{
  fCurrNeuGenCardPair   = 0;

  fCurrDBTable          = 0;
  fCurrDBTableType      = eDbt_Unknown;

  fNeuGenCardPairStack  = new NGCardPairList;
  fNuXSecTableStack     = new DBTableStack<vXSecTableRow>;
  fElDiffXSecTableStack = new DBTableStack<eDiffXSecTableRow>;  
  fSFTableStack         = new DBTableStack<SFTableRow>;
}
//_____________________________________________________________________________
NuVldUserData::~NuVldUserData()
{

}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<vXSecTableRow> * table)
{
  this->DelCurrDBTable();

  fCurrDBTable     = dynamic_cast<DBTableBase *> (table); 
  fCurrDBTableType = eDbt_NuXSec;
}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<eDiffXSecTableRow> * table)
{
  this->DelCurrDBTable();

  fCurrDBTable     = dynamic_cast<DBTableBase *> (table);
  fCurrDBTableType = eDbt_ElDiffXSec;    
}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<SFTableRow> * table)
{
  this->DelCurrDBTable();

  fCurrDBTable     = dynamic_cast<DBTableBase *> (table);
  fCurrDBTableType = eDbt_SF;
}
//_____________________________________________________________________________
void NuVldUserData::DelCurrDBTable(void)
{
  if(fCurrDBTable) delete fCurrDBTable;

  fCurrDBTable = 0;
}
//_____________________________________________________________________________
void NuVldUserData::AddCurrDBTableToStack(string name)
{
  LOG("NuVld", pDEBUG) << "Stacking DBTable<T>: " << name;
  LOG("NuVld", pDEBUG) << "Table-Type: "
                                    << DBTableType::AsString(fCurrDBTableType);

  DBTable<vXSecTableRow> *     nuxsec_table     = 0;
  DBTable<eDiffXSecTableRow> * eldiffxsec_table = 0;
  DBTable<SFTableRow> *        sf_table         = 0;

  if(fCurrDBTable) {
    
    switch(fCurrDBTableType) {
      
    case (eDbt_NuXSec):    

       nuxsec_table = dynamic_cast< DBTable<vXSecTableRow> * > (fCurrDBTable);    
       fNuXSecTableStack->AddDBTable(name, nuxsec_table);

       LOG("NuVld", pDEBUG) << "... done";
       break;
  
    case (eDbt_ElDiffXSec):

       eldiffxsec_table =
                 dynamic_cast< DBTable<eDiffXSecTableRow> * > (fCurrDBTable);
       fElDiffXSecTableStack->AddDBTable(name, eldiffxsec_table);

       LOG("NuVld", pDEBUG) << "... done";
       break;

    case (eDbt_SF):

       sf_table = dynamic_cast< DBTable<SFTableRow> * > (fCurrDBTable);
       fSFTableStack->AddDBTable(name, sf_table);

       LOG("NuVld", pDEBUG) << "... done";
       break;

    default:
       LOG("NuVld", pERROR) << "Could not stack the current DBTable<T>";
    }   
  }
}
//_____________________________________________________________________________
void NuVldUserData::SetStackedDBTableAsCurr(string name)
{
  this->DelCurrDBTable();

  if ( fNuXSecTableStack->Exists(name) ) {
      
       DBTable<vXSecTableRow> * table = this->NuXSecStack()->GetTable(name);
       if(table) fCurrDBTable = new DBTable<vXSecTableRow> (table);

       fCurrDBTableType = eDbt_NuXSec;
  }  

  if ( fElDiffXSecTableStack->Exists(name) ) {

       DBTable<eDiffXSecTableRow> * table =
                                   this->ElDiffXSecStack()->GetTable(name);                                   
       if(table) fCurrDBTable = new DBTable<eDiffXSecTableRow> (table);

       fCurrDBTableType = eDbt_ElDiffXSec;        
  }

  if ( fElDiffXSecTableStack->Exists(name) ) {

       DBTable<SFTableRow> * table = this->SFStack()->GetTable(name);
       if(table) fCurrDBTable = new DBTable<SFTableRow> (table);

       fCurrDBTableType = eDbt_SF;
  }  
}
//_____________________________________________________________________________
void NuVldUserData::DelStackedDBTable(string name)
{
  if (fNuXSecTableStack    ->Exists(name)) fNuXSecTableStack->Erase(name);
  if (fElDiffXSecTableStack->Exists(name)) fElDiffXSecTableStack->Erase(name);
  if (fSFTableStack        ->Exists(name)) fSFTableStack->Erase(name);
}
//_____________________________________________________________________________
DBTable<vXSecTableRow> * NuVldUserData::NuXSec(void)
{
  if(fCurrDBTableType != eDbt_NuXSec) return 0;
  
  DBTable<vXSecTableRow> * table =
                       dynamic_cast< DBTable<vXSecTableRow> * > (fCurrDBTable);  
  return table;
}
//_____________________________________________________________________________
DBTable<eDiffXSecTableRow> * NuVldUserData::ElDiffXSec(void)
{
  if(fCurrDBTableType != eDbt_ElDiffXSec) return 0;

  DBTable<eDiffXSecTableRow> * table =
                  dynamic_cast< DBTable<eDiffXSecTableRow> * > (fCurrDBTable);
  return table;
}
//_____________________________________________________________________________
DBTable<SFTableRow> * NuVldUserData::SF(void)
{
  if(fCurrDBTableType != eDbt_SF) return 0;

  DBTable<SFTableRow> * table =
                         dynamic_cast< DBTable<SFTableRow> * > (fCurrDBTable);
  return table;
}
//_____________________________________________________________________________
unsigned int NuVldUserData::NStackedDBTables(void) const
{
  unsigned int ntables = 0;

  ntables += fNuXSecTableStack     -> GetNTables();
  ntables += fElDiffXSecTableStack -> GetNTables();
  ntables += fSFTableStack         -> GetNTables();

  return ntables;
}
//_____________________________________________________________________________
unsigned int NuVldUserData::NStackedConfigs(void) const
{
  unsigned int nconfigs = 0;

  nconfigs += fNeuGenCardPairStack->NPairs();

  return nconfigs;
}
//_____________________________________________________________________________
void NuVldUserData::ClearStackedConfig(void)
{
  if(fNeuGenCardPairStack) delete fNeuGenCardPairStack;

  fNeuGenCardPairStack = new NGCardPairList;
}
//_____________________________________________________________________________
void NuVldUserData::ClearStackedDBTables(void)
{
  if ( fNuXSecTableStack     ) delete fNuXSecTableStack;
  if ( fElDiffXSecTableStack ) delete fElDiffXSecTableStack;
  if ( fSFTableStack )         delete fSFTableStack;
  
  fNuXSecTableStack     = new DBTableStack<vXSecTableRow>;
  fElDiffXSecTableStack = new DBTableStack<eDiffXSecTableRow>;
  fSFTableStack         = new DBTableStack<SFTableRow>;
}
//_____________________________________________________________________________
const vector<string> * NuVldUserData::GetStackedDBTableNames(void) const
{
  unsigned int itable  = 0;
  unsigned int ntables = this->NStackedDBTables();

  LOG("NuVld", pDEBUG) << "Found " << ntables << " stacked DBTable<T>";

  if(ntables == 0) return 0;

  const vector<string> * names = 0;
  vector<string>::const_iterator iter;
  
  vector<string> * dbtnames = new vector<string>(ntables);

  names = fNuXSecTableStack->GetListOfNames();
     
  for(iter = names->begin(); iter != names->end(); ++iter)
  {
    LOG("NuVld", pDEBUG) << "Adding DBTable<vXSecTableRow> " << *iter;
    (*dbtnames)[itable++] = *iter;
  }
  delete names;                   

  names = fElDiffXSecTableStack->GetListOfNames();

  for(iter = names->begin(); iter != names->end(); ++iter)
  {
    LOG("NuVld", pDEBUG) << "Adding DBTable<eDiffXSecTableRow> " << *iter;
    (*dbtnames)[itable++] = *iter;
  }                   
  delete names;

  names = fSFTableStack->GetListOfNames();

  for(iter = names->begin(); iter != names->end(); ++iter)
  {
    LOG("NuVld", pDEBUG) << "Adding DBTable<SFTableRow> " << *iter;
    (*dbtnames)[itable++] = *iter;
  }
  delete names;

  return dbtnames;
}
//______________________________________________________________________________
void NuVldUserData::SaveStack(string root_file) const
{
  fNuXSecTableStack     -> SaveQueryStringsToFile (root_file, "v_xsec",      "RECREATE");
  fElDiffXSecTableStack -> SaveQueryStringsToFile (root_file, "e_diff_xsec", "UPDATE");
  fSFTableStack         -> SaveQueryStringsToFile (root_file, "sf",          "UPDATE");

  TFile f(root_file.c_str(), "UPDATE");

  fNeuGenCardPairStack -> Write ("generator_config");

  f.Close();
}
//______________________________________________________________________________


