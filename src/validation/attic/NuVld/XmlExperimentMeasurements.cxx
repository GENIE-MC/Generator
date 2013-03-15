//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlExperimentMeasurements.h"

namespace genie {
namespace nuvld {
  
//_______________________________________________________________________________
XmlExperimentMeasurements::XmlExperimentMeasurements()
{
  _data    = new vector<XmlMeasurement *>;
  _spectra = new map<string, XmlBeamFluxSpectrum *>;
}
//_______________________________________________________________________________
XmlExperimentMeasurements::XmlExperimentMeasurements(
                                  const XmlExperimentMeasurements & /*meas_list*/)
{

}
//_______________________________________________________________________________
XmlExperimentMeasurements::~XmlExperimentMeasurements()
{
  delete _data;
  delete _spectra;
}
//_______________________________________________________________________________
void XmlExperimentMeasurements::Add(XmlExperimentInfo * info)
{ 
  _exp_info = info;
}
//_______________________________________________________________________________
void XmlExperimentMeasurements::Add(XmlMeasurement * meas)
{
  _data->push_back(meas);
}
//_______________________________________________________________________________
void XmlExperimentMeasurements::Add(string beam, XmlBeamFluxSpectrum * spectrum)
{
  _spectra->insert(
                 map<string, XmlBeamFluxSpectrum *>::value_type(beam, spectrum) );
}
//_______________________________________________________________________________
const XmlExperimentInfo & XmlExperimentMeasurements::GetXmlExperimentInfo(void) const
{
  return *_exp_info;
}
//_______________________________________________________________________________
const vector<XmlMeasurement *> & XmlExperimentMeasurements::GetXmlMeasurements(void) const
{
  return *_data;
}
//_______________________________________________________________________________
const map<string, XmlBeamFluxSpectrum *> &
                                   XmlExperimentMeasurements::GetSpectra(void) const
{
  return *_spectra;
}
//_______________________________________________________________________________

} // nuvld namespace
} // genie namespace
