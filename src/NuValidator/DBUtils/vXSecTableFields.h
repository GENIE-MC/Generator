//_____________________________________________________________________________
/*!

\class    genie::nuvld::vXSecTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _XSEC_TABLE_FIELDS_H_
#define _XSEC_TABLE_FIELDS_H_

#include "DBUtils/DBTableFields.h"

namespace genie {
namespace nuvld {

class vXSecTableFields : public DBTableFields
{
public:

  vXSecTableFields();
  vXSecTableFields(const DBTableFields * fields);
  virtual ~vXSecTableFields();

ClassDef(vXSecTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif
