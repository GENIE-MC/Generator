//_____________________________________________________________________________
/*!

\class    genie::nuvld::Measurement

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created August 2003
*/
//_____________________________________________________________________________

#ifndef _MEASUREMENT_H_
#define _MEASUREMENT_H_

#include <vector>

#include "MeasurementHeader.h"
#include "Observable.h"
#include "RecordBase.h"

using std::vector;

namespace genie {
namespace nuvld {

class Measurement
{
public:

  Measurement();
  Measurement(const Measurement & meas);
  virtual ~Measurement() { }

  void Add(MeasurementHeader * header);
  void Add(RecordBase * rec);

  const MeasurementHeader &    GetMeasurementHeader(void) const { return *_header; }
  const vector<RecordBase *> & GetDataPoints(void)        const { return *_data;   }

protected:

  Observable_t            _observable;
  MeasurementHeader *     _header;
  vector<RecordBase *> *  _data;
};

} // nuvld namespace
} // genie namespace

#endif // _MEASUREMENT_H_
