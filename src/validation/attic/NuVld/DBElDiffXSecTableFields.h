//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBElDiffXSecTableFields 

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _E_DIFF_XSEC_TABLE_FIELDS_H_
#define _E_DIFF_XSEC_TABLE_FIELDS_H_

#include "ValidationTools/NuVld/DBTableFields.h"

namespace genie {
namespace nuvld {

class DBElDiffXSecTableFields : public DBTableFields
{
public:

  DBElDiffXSecTableFields();
  DBElDiffXSecTableFields(const DBTableFields * fields);
  virtual ~DBElDiffXSecTableFields();

ClassDef(DBElDiffXSecTableFields, 1)
};

} // nuvld namespace
} // genie namespace

#endif // _E_DIFF_XSEC_TABLE_FIELDS_H_
