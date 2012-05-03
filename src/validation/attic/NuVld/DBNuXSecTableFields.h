//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBNuXSecTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _XSEC_TABLE_FIELDS_H_
#define _XSEC_TABLE_FIELDS_H_

#include  "ValidationTools/NuVld/DBTableFields.h"

namespace genie {
namespace nuvld {

class DBNuXSecTableFields : public DBTableFields
{
public:

  DBNuXSecTableFields();
  DBNuXSecTableFields(const DBTableFields * fields);
  virtual ~DBNuXSecTableFields();

ClassDef(DBNuXSecTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif
