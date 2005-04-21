//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBI

\brief    Data Base Interface

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#ifndef _NUVLD_KERNEL_H_
#define _NUVLD_KERNEL_H_

#include <TSQLServer.h>

#include "DBUtils/DBStatus.h"
#include "DBUtils/DBTable.h"
#include "DBUtils/DBQueryString.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "DBUtils/SFTableRow.h"
#include "XmlParser/XmlDataSet.h"

namespace genie {
namespace nuvld {

class DBI {

public:

  DBI();
  DBI(TSQLServer * sql_server);
  ~DBI();

  DBStatus_t  UploadXML (const XmlDataSet & data);
  DBStatus_t  FillTable (genie::nuvld::DBTable<vXSecTableRow> * t,     const DBQueryString & query_string);
  DBStatus_t  FillTable (genie::nuvld::DBTable<eDiffXSecTableRow> * t, const DBQueryString & query_string);
  DBStatus_t  FillTable (genie::nuvld::DBTable<SFTableRow> * t,        const DBQueryString & query_string);
  
private:

  bool HaveConnection(void) const;
  void AddMeasurementIdList(DBTableBase * table);

  TSQLServer * fSqlServer;
};

} // nuvld namespace
} // genie namespace

#endif

