//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _DB_TABLE_ROW_H_
#define _DB_TABLE_ROW_H_

#include <vector>
#include <string>
#include <ostream>

#include <TObject.h>
#include <TSQLRow.h>

#include "DBUtils/DBTableFields.h"

using std::vector;
using std::string;
using std::ostream;

namespace genie {
namespace nuvld {

class DBTableRow : public TObject
{
public:

  string Field (unsigned int field_pos) const;
  string Field (string field_name) const;

  void Print(ostream & stream) const;

protected:

  DBTableRow();
  DBTableRow(const DBTableFields * fields, TSQLRow * row);
  DBTableRow(const DBTableRow * row);
  ~DBTableRow();

  vector<string>        _row;
  const DBTableFields * _fields;

ClassDef(DBTableRow, 1)
};

} // nuvld namespace
} // genie namespace

#endif
