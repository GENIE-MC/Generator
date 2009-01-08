//_____________________________________________________________________________
/*!

\class    genie::nuvld::ExperimentMeasurements

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "ExperimentMeasurements.h"

namespace genie {
namespace nuvld {
  
//_______________________________________________________________________________
ExperimentMeasurements::ExperimentMeasurements()
{
  _data    = new vector<Measurement *>;
  _spectra = new map<string, BeamFluxSpectrum *>;
}
//_______________________________________________________________________________
ExperimentMeasurements::ExperimentMeasurements(
                                  const ExperimentMeasurements & /*meas_list*/)
{

}
//_______________________________________________________________________________
ExperimentMeasurements::~ExperimentMeasurements()
{
  delete _data;
  delete _spectra;
}
//_______________________________________________________________________________
void ExperimentMeasurements::Add(ExperimentInfo * info)
{ 
  _exp_info = info;
}
//_______________________________________________________________________________
void ExperimentMeasurements::Add(Measurement * meas)
{
  _data->push_back(meas);
}
//_______________________________________________________________________________
void ExperimentMeasurements::Add(string beam, BeamFluxSpectrum * spectrum)
{
  _spectra->insert(
                 map<string, BeamFluxSpectrum *>::value_type(beam, spectrum) );
}
//_______________________________________________________________________________
const ExperimentInfo & ExperimentMeasurements::GetExperimentInfo(void) const
{
  return *_exp_info;
}
//_______________________________________________________________________________
const vector<Measurement *> & ExperimentMeasurements::GetMeasurements(void) const
{
  return *_data;
}
//_______________________________________________________________________________
const map<string, BeamFluxSpectrum *> &
                                   ExperimentMeasurements::GetSpectra(void) const
{
  return *_spectra;
}
//_______________________________________________________________________________

} // nuvld namespace
} // genie namespace
