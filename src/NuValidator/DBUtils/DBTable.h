//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTable

\brief    NuVld DBase SQL Query Result

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_H_
#define _DB_TABLE_H_

#include <iostream>
#include <vector>
#include <string>

#include <TDirectory.h>
#include <TGraphAsymmErrors.h>

#include "DBUtils/DBTableBase.h"
#include "DBUtils/MeasurementIdList.h"
#include "DBUtils/MultiGraph.h"

using std::ostream;
using std::vector;
using std::string;

namespace genie {
namespace nuvld {

template<class T> class DBTable;
//template<class T> ostream & operator << (ostream & stream, const DBTable<T> & table);

template<class T> class DBTable : public DBTableBase
{
public:

   friend class DBI;
   //friend ostream & operator << <T>(ostream & stream, const DBTable<T> & table);
   
   DBTable();
   DBTable(const genie::nuvld::DBTable<T> * table);
   ~DBTable();
   
   int   NRows            (void)  const;
   void  MergeWithTable   (DBTableBase * table);
   
   const DBTableRow *        Row     (int irow) const;
   const DBTableFields *     Fields      (void) const { return _fields;      }
   const vector<T *> &       Rows        (void) const { return _table;       }
   const MeasurementIdList * IdList      (void) const { return _id_list;     }
   const DBQueryString *     QueryString (void) const { return _query_string;}
   
   DBTable * Subset(string experiment, string measurement_tag) const;

   TGraphAsymmErrors * GetGraph      (const char * opt, const char * var = 0);
   MultiGraph *        GetMultiGraph (const char * opt, const char * var = 0);

   void SaveQueryStringToFile (TDirectory * dir, string name) const;
   //void LoadFromFile (TDirectory * dir, const char * name);

private:

   void AddRow               (DBTableRow * row);
   void SetQueryString       (const DBQueryString & query_string);
   void SetMeasurementIdList (MeasurementIdList * id_list);

   T *                   _db_table_row;
   DBTableFields *       _fields;
   vector<T *>           _table;
   MeasurementIdList *   _id_list;
   DBQueryString *       _query_string;
};

} // nuvld namespace
} // genie namespace

#endif
