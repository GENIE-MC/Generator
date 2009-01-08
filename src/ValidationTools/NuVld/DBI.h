//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBI

\brief    Data Base Interface

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#ifndef _DBI_H_
#define _DBI_H_

#include <TSQLServer.h>

#include "ValidationTools/NuVld/DBStatus.h"
#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBQueryString.h"
#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBSFTableRow.h"
#include "ValidationTools/NuVld/XmlDataSet.h"

namespace genie {
namespace nuvld {

class DBI {

public:
  DBI();
  DBI(TSQLServer * sql_server);
  ~DBI();

  DBStatus_t  UploadXML (const XmlDataSet & data);
  
  DBStatus_t  FillTable (genie::nuvld::DBTable<DBNuXSecTableRow> * t,     const DBQueryString & q);
  DBStatus_t  FillTable (genie::nuvld::DBTable<DBElDiffXSecTableRow> * t, const DBQueryString & q);
  DBStatus_t  FillTable (genie::nuvld::DBTable<DBSFTableRow> * t,         const DBQueryString & q);
  
private:
  bool HaveConnection(void) const;
  void AddDBMeasurementIdList(DBTableBase * table);

  TSQLServer * fSqlServer;
};

} // nuvld namespace
} // genie namespace

#endif

