//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBMeasurementId

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _MEASUREMENT_ID_H_
#define _MEASUREMENT_ID_H_

#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::vector;
using std::ostream;

namespace genie {
namespace nuvld {

class DBMeasurementId 
{
public:

  friend class DBI;

  friend ostream & operator << (ostream & stream, const DBMeasurementId & id);

  DBMeasurementId();
  DBMeasurementId(const DBMeasurementId * mid);
  DBMeasurementId(string experiment, string measurement_tag);
  ~DBMeasurementId();

  unsigned int NRefs           (void)                  const;
  string       Experiment      (void)                  const;
  string       XmlMeasurementTag  (void)                  const;
  string       Reference       (unsigned int iref = 0) const;
  string       Author          (unsigned int iref = 0) const;
  string       Journal         (unsigned int iref = 0) const;
  string       Year            (unsigned int iref = 0) const;
  void         Print           (ostream & stream)      const;

private:

  string         _experiment;
  string         _measurement_tag;
  vector<string> _author;
  vector<string> _journal;
  vector<string> _year;
};

} // nuvld namespace
} // genie namespace

#endif
