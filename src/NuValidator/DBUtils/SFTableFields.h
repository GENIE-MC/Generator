//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFTableFields

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _SF_TABLE_FIELDS_H_
#define _SF_TABLE_FIELDS_H_

#include "DBUtils/DBTableFields.h"

namespace genie {
namespace nuvld {

class SFTableFields : public DBTableFields
{
public:

  SFTableFields();
  SFTableFields(const DBTableFields * fields);
  virtual ~SFTableFields();

ClassDef(SFTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif
