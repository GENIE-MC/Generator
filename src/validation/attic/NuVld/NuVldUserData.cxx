//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 21, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TFile.h>

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/NuVldUserData.h"

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
//fCurrNeuGenCardPair   = 0;

  fCurrDBTable          = 0;
  fCurrDBTableType      = eDbt_Unknown;

  fNuXSecTableStack     = new DBTableStack<DBNuXSecTableRow>;
  fElDiffXSecTableStack = new DBTableStack<DBElDiffXSecTableRow>;  
  fSFTableStack         = new DBTableStack<DBSFTableRow>;
}
//_____________________________________________________________________________
NuVldUserData::~NuVldUserData()
{

}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<DBNuXSecTableRow> * table)
{
  this->DelCurrDBTable();

  fCurrDBTable     = dynamic_cast<DBTableBase *> (table); 
  fCurrDBTableType = eDbt_NuXSec;
}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<DBElDiffXSecTableRow> * table)
{
  this->DelCurrDBTable();

  fCurrDBTable     = dynamic_cast<DBTableBase *> (table);
  fCurrDBTableType = eDbt_ElDiffXSec;    
}
//_____________________________________________________________________________
void NuVldUserData::SetCurrDBTable(DBTable<DBSFTableRow> * table)
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

  DBTable<DBNuXSecTableRow> *     nuxsec_table     = 0;
  DBTable<DBElDiffXSecTableRow> * eldiffxsec_table = 0;
  DBTable<DBSFTableRow> *         sf_table         = 0;

  if(fCurrDBTable) {
    
    switch(fCurrDBTableType) {
      
    case (eDbt_NuXSec):    

       nuxsec_table = dynamic_cast< DBTable<DBNuXSecTableRow> * > (fCurrDBTable);    
       fNuXSecTableStack->AddDBTable(name, nuxsec_table);

       LOG("NuVld", pDEBUG) << "... done";
       break;
  
    case (eDbt_ElDiffXSec):

       eldiffxsec_table =
                 dynamic_cast< DBTable<DBElDiffXSecTableRow> * > (fCurrDBTable);
       fElDiffXSecTableStack->AddDBTable(name, eldiffxsec_table);

       LOG("NuVld", pDEBUG) << "... done";
       break;

    case (eDbt_SF):

       sf_table = dynamic_cast< DBTable<DBSFTableRow> * > (fCurrDBTable);
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
       DBTable<DBNuXSecTableRow> * table = this->NuXSecStack()->GetTable(name);
       if(table) fCurrDBTable = new DBTable<DBNuXSecTableRow> (table);
       fCurrDBTableType = eDbt_NuXSec;
  }  

  if ( fElDiffXSecTableStack->Exists(name) ) {
       DBTable<DBElDiffXSecTableRow> * table = this->ElDiffXSecStack()->GetTable(name);                                   
       if(table) fCurrDBTable = new DBTable<DBElDiffXSecTableRow> (table);
       fCurrDBTableType = eDbt_ElDiffXSec;        
  }

  if ( fElDiffXSecTableStack->Exists(name) ) {
       DBTable<DBSFTableRow> * table = this->SFStack()->GetTable(name);
       if(table) fCurrDBTable = new DBTable<DBSFTableRow> (table);
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
DBTable<DBNuXSecTableRow> * NuVldUserData::NuXSec(void)
{
  if(fCurrDBTableType != eDbt_NuXSec) return 0;
  
  DBTable<DBNuXSecTableRow> * table =
       dynamic_cast< DBTable<DBNuXSecTableRow> * > (fCurrDBTable);  
  return table;
}
//_____________________________________________________________________________
DBTable<DBElDiffXSecTableRow> * NuVldUserData::ElDiffXSec(void)
{
  if(fCurrDBTableType != eDbt_ElDiffXSec) return 0;

  DBTable<DBElDiffXSecTableRow> * table =
      dynamic_cast< DBTable<DBElDiffXSecTableRow> * > (fCurrDBTable);
  return table;
}
//_____________________________________________________________________________
DBTable<DBSFTableRow> * NuVldUserData::SF(void)
{
  if(fCurrDBTableType != eDbt_SF) return 0;

  DBTable<DBSFTableRow> * table =
     dynamic_cast< DBTable<DBSFTableRow> * > (fCurrDBTable);
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

//nconfigs += fNeuGenCardPairStack->NPairs();

  return nconfigs;
}
//_____________________________________________________________________________
void NuVldUserData::ClearStackedConfig(void)
{
/*
  if(fNeuGenCardPairStack) delete fNeuGenCardPairStack;

  fNeuGenCardPairStack = new NGCardPairList;
*/
}
//_____________________________________________________________________________
void NuVldUserData::ClearStackedDBTables(void)
{
  if ( fNuXSecTableStack     ) delete fNuXSecTableStack;
  if ( fElDiffXSecTableStack ) delete fElDiffXSecTableStack;
  if ( fSFTableStack )         delete fSFTableStack;
  
  fNuXSecTableStack     = new DBTableStack<DBNuXSecTableRow>;
  fElDiffXSecTableStack = new DBTableStack<DBElDiffXSecTableRow>;
  fSFTableStack         = new DBTableStack<DBSFTableRow>;
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
    LOG("NuVld", pDEBUG) << "Adding DBTable<DBNuXSecTableRow> " << *iter;
    (*dbtnames)[itable++] = *iter;
  }
  delete names;                   

  names = fElDiffXSecTableStack->GetListOfNames();

  for(iter = names->begin(); iter != names->end(); ++iter)
  {
    LOG("NuVld", pDEBUG) << "Adding DBTable<DBElDiffXSecTableRow> " << *iter;
    (*dbtnames)[itable++] = *iter;
  }                   
  delete names;

  names = fSFTableStack->GetListOfNames();

  for(iter = names->begin(); iter != names->end(); ++iter)
  {
    LOG("NuVld", pDEBUG) << "Adding DBTable<DBSFTableRow> " << *iter;
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

//fNeuGenCardPairStack -> Write ("generator_config");

  f.Close();
}
//______________________________________________________________________________


