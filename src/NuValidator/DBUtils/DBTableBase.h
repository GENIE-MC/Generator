//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableBase

\brief    

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_BASE_H_
#define _DB_TABLE_BASE_H_

#include <vector>
#include <string>

#include "DBUtils/DBTableFields.h"
#include "DBUtils/DBTableRow.h"
#include "DBUtils/DBQueryString.h"
#include "DBUtils/MeasurementIdList.h"

using std::vector;
using std::string;

namespace genie {
namespace nuvld {

class DBTableBase
{
public:

   DBTableBase() { }
   virtual ~DBTableBase() { };
   
   virtual int   NRows          (void)  const              = 0;
   virtual void  AddRow         (DBTableRow * row)         = 0;
   virtual void  MergeWithTable (DBTableBase * table)      = 0;
   virtual void  SetQueryString (const DBQueryString & qs) = 0;

   
   virtual const DBTableRow *        Row     (int irow) const = 0;
   virtual const DBTableFields *     Fields      (void) const = 0;
   virtual const MeasurementIdList * IdList      (void) const = 0;
   virtual const DBQueryString *     QueryString (void) const = 0;
   
//   virtual DBTableBase * Subset(string exp, string mtag) const = 0;

   virtual void SetMeasurementIdList(MeasurementIdList * id_list) = 0;
};

} // nuvld namespace
} // genie namespace

#endif
