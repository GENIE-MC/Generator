//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDiffXSecTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _E_DIFF_XSEC_TABLE_FIELDS_H_
#define _E_DIFF_XSEC_TABLE_FIELDS_H_

#include "DBUtils/DBTableFields.h"

namespace genie {
namespace nuvld {

class eDiffXSecTableFields : public DBTableFields
{
public:

  eDiffXSecTableFields();
  eDiffXSecTableFields(const DBTableFields * fields);
  virtual ~eDiffXSecTableFields();

ClassDef(eDiffXSecTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif // _E_DIFF_XSEC_TABLE_FIELDS_H_
