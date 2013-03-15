//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlMeasurement

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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
