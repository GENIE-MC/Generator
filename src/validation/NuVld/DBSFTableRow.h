//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBSFTableRow 

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _SF_TABLE_ROW_H_
#define _SF_TABLE_ROW_H_

#include <string>
#include <ostream>

#include <TSQLRow.h>

#include  "ValidationTools/NuVld/DBTableRow.h"

using std::string;
using std::ostream;

namespace genie {
namespace nuvld {

class DBSFTableRow : public DBTableRow
{
public:

  friend class DBI;

  friend ostream & operator << (ostream & stream, const DBSFTableRow & row);

  string Experiment     (void) const;  
  string XmlMeasurementTag (void) const;
  double SF             (void) const; ///< structure function
  string R              (void) const; ///< longitud./transv. structure function
  string p              (void) const; ///< Proton/Neutron/Nuclear
  double StatErrP       (void) const;
  double StatErrM       (void) const;
  double SystErrP       (void) const;
  double SystErrM       (void) const;
  double ErrP           (void) const;
  double ErrM           (void) const;
  double x              (void) const;
  double Q2             (void) const;

  void   Print (ostream & stream) const;

private:

  DBSFTableRow ();
  DBSFTableRow (const DBTableRow * db_row);
  DBSFTableRow (TSQLRow * row);
  virtual ~DBSFTableRow ();

ClassDef(DBSFTableRow, 1)
};

} // nuvld namespace
} // genie namespace

#endif
