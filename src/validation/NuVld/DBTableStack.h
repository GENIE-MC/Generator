//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableStack

\brief    A collection of named DBTable<T> objects

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_STACK_H_
#define _DB_TABLE_STACK_H_

#include <map>
#include <vector>
#include <string>

#include  "ValidationTools/NuVld/DBTable.h"

using std::map;
using std::vector;
using std::string;

namespace genie {
namespace nuvld {

template<class T> class DBTableStack;

//template<class T> ostream & operator
//                      << (ostream & stream, const DBTableStack<T> & mtable);

template<class T> class DBTableStack 
{
public:

//  friend ostream & operator
//                   << <T>(ostream & stream, const DBTableStack<T> & mtable);

  DBTableStack();
  DBTableStack(const DBTableStack * mt);
  ~DBTableStack();

  void                   AddDBTable     (string name, genie::nuvld::DBTable<T> * table);
  void                   Merge          (const DBTableStack * mt);
  void                   Erase          (string name);
  const vector<string> * GetListOfNames (void)        const;
  DBTable<T> *           GetTable       (string name) const;
  unsigned int           GetNTables     (void)        const;
  bool                   Exists         (string name) const;

  void SaveQueryStringsToFile (string filename, string name, string option) const;

private:

  map<string, genie::nuvld::DBTable<T> * > _db_table_map;
};


} // nuvld namespace
} // genie namespace

#endif
