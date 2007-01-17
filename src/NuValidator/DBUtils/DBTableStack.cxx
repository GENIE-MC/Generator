//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableStack

\brief    A collection of named DBTable<T> objects

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#include <iostream>

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>

#include "DBUtils/DBTableStack.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "DBUtils/SFTableRow.h"
#include "Messenger/Messenger.h"

using std::endl;
using std::cout;
using std::ostringstream;

namespace genie {
namespace nuvld {

//template ostream & operator
//           << (ostream & stream, const DBTableStack<vXSecTableRow> & table);
//template ostream & operator
//       << (ostream & stream, const DBTableStack<eDiffXSecTableRow> & table);
//template ostream & operator
//              << (ostream & stream, const DBTableStack<SFTableRow> & table);

//____________________________________________________________________________
/*
template<class T>
    ostream & operator << (ostream & stream, const DBTableStack<T> & mtable)
{
  typename map<string, DBTable<T> *>::const_iterator iter;

  stream << "Contained xsec_tables: " << endl;
  for(iter = mtable._db_table_map.begin();
                                 iter != mtable._db_table_map.end(); ++iter) {

    stream << "xsec_table: " << iter->first << " with "
                                  << iter->second->NRows() << " rows" << endl;
  }
  return stream;
}*/
//______________________________________________________________________________
template<class T> DBTableStack<T>::DBTableStack()
{

}
//______________________________________________________________________________
template<class T> DBTableStack<T>::DBTableStack(const DBTableStack<T> * mt)
{
  typename map<string, DBTable<T> *>::const_iterator mt_iter;

  for(mt_iter = mt->_db_table_map.begin();
                               mt_iter != mt->_db_table_map.end(); ++mt_iter) {

     string       name  =  mt_iter->first;
     DBTable<T> * table =  mt_iter->second;

     this->AddDBTable(name, table);
  }
}
//______________________________________________________________________________
template<class T> DBTableStack<T>::~DBTableStack()
{

}
//______________________________________________________________________________
template<class T> void DBTableStack<T>::AddDBTable(
                                               string name, DBTable<T> * table)
{
  DBTable<T> * cloned_table = new DBTable<T>(table);

  _db_table_map.insert(
            typename map<string, DBTable<T> *>::value_type(name,cloned_table) );
}
//______________________________________________________________________________
template<class T> void DBTableStack<T>::Merge(const DBTableStack * mt)
{
  typename map<string, DBTable<T> *>::const_iterator mt_iter;

  for(mt_iter = mt->_db_table_map.begin();
                               mt_iter != mt->_db_table_map.end(); ++mt_iter) {

     string        name =  mt_iter->first;
     DBTable<T> * table =  mt_iter->second;

     this->AddDBTable(name, table);
  }
}
//______________________________________________________________________________
template<class T> void DBTableStack<T>::Erase(string name)
{
   if( this->Exists(name) ) _db_table_map.erase(name);
}
//______________________________________________________________________________
template<class T>
              const vector<string> * DBTableStack<T>::GetListOfNames(void) const
{
   typename map<string, DBTable<T> *>::const_iterator iter;

   vector<string> * names = new vector<string>;

   for(iter = _db_table_map.begin();
                            iter != _db_table_map.end(); ++iter)
                                                names->push_back( iter->first );
   return names;
}
//______________________________________________________________________________
template<class T> DBTable<T> * DBTableStack<T>::GetTable(string name) const
{
  if( _db_table_map.count(name) == 1 ) {

   typename map<string, DBTable<T> *>::const_iterator iter =
                                                       _db_table_map.find(name);
   return iter->second;

  } else return 0;
}
//______________________________________________________________________________
template<class T> unsigned int DBTableStack<T>::GetNTables(void) const
{
  return _db_table_map.size();
}
//______________________________________________________________________________
template<class T> bool DBTableStack<T>::Exists(string name) const
{
  return (_db_table_map.count(name) == 1);
}
//______________________________________________________________________________
template<class T> void DBTableStack<T>::SaveQueryStringsToFile(
                              string filename, string name, string option) const
{
  typename map<string, DBTable<T> *>::const_iterator iter;

  TFile f(filename.c_str(), option.c_str());

  TDirectory * top_dir = f.mkdir(name.c_str());

  top_dir->cd();

  for(iter = _db_table_map.begin(); iter != _db_table_map.end(); ++iter) {

    LOG("NuVld", pDEBUG)
                << "Saving the query that creates the table: " << iter->first;

    iter->second->SaveQueryStringToFile(top_dir, iter->first);
  }
}
//____________________________________________________________________________

// template specializations:
//
template class DBTableStack<vXSecTableRow>;
template class DBTableStack<eDiffXSecTableRow>;
template class DBTableStack<SFTableRow>;

//____________________________________________________________________________

} // nuvld namespace
} // genie namespace

