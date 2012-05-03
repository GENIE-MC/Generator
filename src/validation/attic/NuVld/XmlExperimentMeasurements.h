//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlExperimentMeasurements

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _EXPERIMENT_MEASUREMENTS_H_
#define _EXPERIMENT_MEASUREMENTS_H_

#include <map>
#include <vector>

#include "ValidationTools/NuVld/XmlExperimentInfo.h"
#include "ValidationTools/NuVld/XmlBeamFluxSpectrum.h"
#include "ValidationTools/NuVld/XmlMeasurement.h"

using std::map;
using std::vector;

namespace genie {
namespace nuvld {
  
class XmlExperimentMeasurements
{
public:

  XmlExperimentMeasurements();
  XmlExperimentMeasurements(const XmlExperimentMeasurements & meas_list);
  ~XmlExperimentMeasurements();

  void Add (XmlMeasurement * meas);
  void Add (XmlExperimentInfo * info);
  void Add (string beam, XmlBeamFluxSpectrum * spectrum);

  const XmlExperimentInfo &  GetXmlExperimentInfo(void) const;
  const vector<XmlMeasurement *> & GetXmlMeasurements(void) const;    
  const map<string, XmlBeamFluxSpectrum *> & GetSpectra(void) const;

private:

  XmlExperimentInfo *                   _exp_info;
  vector<XmlMeasurement *> *            _data;    
  map<string, XmlBeamFluxSpectrum *> *  _spectra;    
};

} // nuvld namespace
} // genie namespace

#endif // _EXPERIMENT_MEASUREMENTS_H_
