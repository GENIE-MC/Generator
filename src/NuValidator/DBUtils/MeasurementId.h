//_____________________________________________________________________________
/*!

\class    genie::nuvld::MeasurementId

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

class MeasurementId 
{
public:

  friend class DBI;

  friend ostream & operator << (ostream & stream, const MeasurementId & id);

  MeasurementId();
  MeasurementId(const MeasurementId * mid);
  MeasurementId(string experiment, string measurement_tag);
  ~MeasurementId();

  unsigned int NRefs           (void)                  const;
  string       Experiment      (void)                  const;
  string       MeasurementTag  (void)                  const;
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
