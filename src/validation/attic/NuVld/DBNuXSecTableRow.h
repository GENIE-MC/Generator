//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBNuXSecTableRow

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _XSEC_TABLE_ROW_H_
#define _XSEC_TABLE_ROW_H_

#include <string>
#include <ostream>

#include <TSQLRow.h>

#include  "ValidationTools/NuVld/DBTableRow.h"

using std::string;
using std::ostream;

namespace genie {
namespace nuvld {

class DBNuXSecTableRow : public DBTableRow
{
public:

  friend class DBI;

  friend ostream & operator << (ostream & stream, const DBNuXSecTableRow & row);

  string Experiment     (void) const;
  string XmlMeasurementTag (void) const;
  string XSecNorm       (void) const;
  string XSecUnits      (void) const;
  string StatErrType    (void) const;
  string SystErrType    (void) const;
  double XSec           (void) const;
  double StatErrP       (void) const;
  double StatErrM       (void) const;
  double SystErrP       (void) const;
  double SystErrM       (void) const;
  double ErrP           (void) const;
  double ErrM           (void) const;
  double E              (void) const;
  double Emin           (void) const;
  double Emax           (void) const;
  double dE             (void) const;

  void Print (ostream & stream) const;

private:

  DBNuXSecTableRow();
  DBNuXSecTableRow(const DBTableRow * db_row);
  DBNuXSecTableRow(TSQLRow * row);
  virtual ~DBNuXSecTableRow();

  void   UndoXSecNorm     (double & xsec, double E, string norm) const;
  void   ApplyUnitsFactor (double & xsec, string factor) const;
  bool   IsPercentageErr  (string err_type) const;
ClassDef(DBNuXSecTableRow, 1)
};

} // nuvld namespace
} // genie namespace

#endif
