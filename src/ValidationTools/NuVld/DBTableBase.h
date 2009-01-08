//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableBase

\brief    

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_BASE_H_
#define _DB_TABLE_BASE_H_

#include <vector>
#include <string>

#include  "ValidationTools/NuVld/DBTableFields.h"
#include  "ValidationTools/NuVld/DBTableRow.h"
#include  "ValidationTools/NuVld/DBQueryString.h"
#include  "ValidationTools/NuVld/DBMeasurementIdList.h"

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

   virtual const DBTableRow *          Row     (int irow) const = 0;
   virtual const DBTableFields   *     Fields      (void) const = 0;
   virtual const DBMeasurementIdList * IdList      (void) const = 0;
   virtual const DBQueryString *       QueryString (void) const = 0;
   
//   virtual DBTableBase * Subset(string exp, string mtag) const = 0;

   virtual void SetDBMeasurementIdList(DBMeasurementIdList * id_list) = 0;
};

} // nuvld namespace
} // genie namespace

#endif
