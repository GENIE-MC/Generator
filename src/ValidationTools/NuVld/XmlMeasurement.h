//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlMeasurement

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created August 2003
*/
//_____________________________________________________________________________

#ifndef _MEASUREMENT_H_
#define _MEASUREMENT_H_

#include <vector>

#include "ValidationTools/NuVld/XmlMeasurementHeader.h"
#include "ValidationTools/NuVld/XmlObservable.h"
#include "ValidationTools/NuVld/XmlRecordBase.h"

using std::vector;

namespace genie {
namespace nuvld {

class XmlMeasurement
{
public:

  XmlMeasurement();
  XmlMeasurement(const XmlMeasurement & meas);
  virtual ~XmlMeasurement() { }

  void Add(XmlMeasurementHeader * header);
  void Add(XmlRecordBase * rec);

  const XmlMeasurementHeader &    GetXmlMeasurementHeader(void) const { return *_header; }
  const vector<XmlRecordBase *> & GetDataPoints(void)        const { return *_data;   }

protected:

  XmlObservable_t            _observable;
  XmlMeasurementHeader *     _header;
  vector<XmlRecordBase *> *  _data;
};

} // nuvld namespace
} // genie namespace

#endif // _MEASUREMENT_H_
