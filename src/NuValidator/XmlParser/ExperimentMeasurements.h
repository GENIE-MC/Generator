//_____________________________________________________________________________
/*!

\class    genie::nuvld::ExperimentMeasurements

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _EXPERIMENT_MEASUREMENTS_H_
#define _EXPERIMENT_MEASUREMENTS_H_

#include <map>
#include <vector>

#include "ExperimentInfo.h"
#include "BeamFluxSpectrum.h"
#include "Measurement.h"

using std::map;
using std::vector;

namespace genie {
namespace nuvld {
  
class ExperimentMeasurements
{
public:

  ExperimentMeasurements();
  ExperimentMeasurements(const ExperimentMeasurements & meas_list);
  ~ExperimentMeasurements();

  void Add (Measurement * meas);
  void Add (ExperimentInfo * info);
  void Add (string beam, BeamFluxSpectrum * spectrum);

  const ExperimentInfo &  GetExperimentInfo(void) const;
  const vector<Measurement *> & GetMeasurements(void) const;    
  const map<string, BeamFluxSpectrum *> & GetSpectra(void) const;

private:

  ExperimentInfo *                   _exp_info;
  vector<Measurement *> *            _data;    
  map<string, BeamFluxSpectrum *> *  _spectra;    
};

} // nuvld namespace
} // genie namespace

#endif // _EXPERIMENT_MEASUREMENTS_H_
