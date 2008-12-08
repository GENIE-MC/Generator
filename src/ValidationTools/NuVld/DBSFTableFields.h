//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBSFTableFields "ValidationTools/NuVld/

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _SF_TABLE_FIELDS_H_
#define _SF_TABLE_FIELDS_H_

#include  "ValidationTools/NuVld/DBTableFields.h"

namespace genie {
namespace nuvld {

class DBSFTableFields: public DBTableFields
{
public:

  DBSFTableFields();
  DBSFTableFields(const DBTableFields * fields);
  virtual ~DBSFTableFields();

ClassDef(DBSFTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif
