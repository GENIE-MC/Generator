//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <sstream>
#include <string>
#include <map>

#include "ValidationTools/NuVld/DBXmlUploader.h"
#include "Messenger/Messenger.h"

using namespace genie::nuvld;

//_______________________________________________________________________________
DBXmlUploader::DBXmlUploader(TSQLServer * server)
{
 _sql_server = server;
}
//_______________________________________________________________________________
DBXmlUploader::~DBXmlUploader()
{

}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::Upload(const XmlDataSet & nuscat_data) const
{
 const map<string, XmlExperimentMeasurements *> & all_data = nuscat_data.Get();

 map<string, XmlExperimentMeasurements *>::const_iterator exp_iter;

 for(exp_iter = all_data.begin(); exp_iter != all_data.end(); ++exp_iter) {

   const XmlExperimentInfo & exp_info = exp_iter->second->GetXmlExperimentInfo();
   const map<string, XmlBeamFluxSpectrum *> & spectra = exp_iter->second->GetSpectra();

   this->UploadExpInfo( exp_info );
   this->UploadFluxSpectra( exp_info, spectra);

   const vector<XmlMeasurement *> & meas = exp_iter->second->GetXmlMeasurements();

   vector<XmlMeasurement *>::const_iterator meas_iter;

   for(meas_iter = meas.begin(); meas_iter != meas.end(); ++meas_iter) {

       const XmlMeasurementHeader & m_header = (*meas_iter)->GetXmlMeasurementHeader();

       this->UploadReferences( exp_info, m_header );
       this->UploadMeasHdr   ( exp_info, m_header );

       if( m_header.Observable().find("TOT_XSEC") == 0  ||
           m_header.Observable().find("QES_XSEC") == 0  ||
           m_header.Observable().find("SPP_XSEC") == 0  ||
           m_header.Observable().find("COH_XSEC") == 0  ||
           m_header.Observable().find("MPP_XSEC") == 0  ) 
       {
          this->UploadNuXSec(exp_info, *meas_iter);
       }
       else
       if( m_header.Observable().find("ELEC_PXSEC") == 0 )
       {
          this->UploadElDiffXSec(exp_info, *meas_iter);
       }
       else
       if( m_header.Observable().find("F2") == 0 ||
           m_header.Observable().find("xF3") == 0 )
       {
          this->UploadSF(exp_info, *meas_iter);
       }
       else {
          SLOG("NuVld", pERROR) << "Unrecognized observable!";
       }
   }
 }
 return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadExpInfo(const XmlExperimentInfo & exp_info) const
{
 ostringstream sql_string;

 sql_string << "INSERT INTO EXP_INFO VALUES "
            << "(\""  << exp_info.Name()           << "\""
            << ",\""  << exp_info.Comment()        << "\""
            << ",\""  << exp_info.Facility()       << "\""
            << ",\""  << exp_info.Detector()       << "\""
            << ",\""  << exp_info.Beam()           << "\""
            << ",\""  << exp_info.Target()         << "\""
            << ","    << exp_info.YearBegin()
            << ","    << exp_info.YearEnd()
            << ",\""  << exp_info.Exposure()       << "\""
            << ",\""  << exp_info.ExposureUnits()  << "\""
            << ","    << exp_info.EnergyMin()
            << ","    << exp_info.EnergyMax()
            << ",\""  << exp_info.EnergyUnits()    << "\""
            << ",\""  << exp_info.EnergyFrame()    << "\""
            << ");";

  SLOG("NuVld", pINFO)
           << "SQL string to be sent to the DBase: " << sql_string.str().c_str();

  _sql_server->Query( sql_string.str().c_str() );

  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadFluxSpectra(
                                         const XmlExperimentInfo & exp_info,
                         const map<string, XmlBeamFluxSpectrum *> & spectra) const
{
  //-- loop over all flux spectra for this experiment

  map<string, XmlBeamFluxSpectrum *>::const_iterator spectrum_iter;
  for(spectrum_iter = spectra.begin();
                               spectrum_iter != spectra.end(); ++spectrum_iter) {

       const string beam = spectrum_iter->first;
       const vector<XmlBeamFluxBin *> & spectrum_bins
                                          = spectrum_iter->second->GetSpectrum();

       vector<XmlBeamFluxBin *>::const_iterator bin;

       //-- loop over all energy bins for the current spectrum

       for(bin = spectrum_bins.begin(); bin != spectrum_bins.end(); ++bin) {

              ostringstream sql_string;

              sql_string
                  << "INSERT INTO BEAM_FLUX VALUES "
                  << "(\""  << exp_info.Name()        << "\""
                  << ",\""  << beam                   << "\""
                  << ","    << (*bin)->MeanEnergy()
                  << ","    << (*bin)->MinEnergy()
                  << ","    << (*bin)->MaxEnergy()
                  << ",\""  << spectrum_iter->second->EnergyUnits() << "\""
                  << ",\""  << spectrum_iter->second->EnergyFrame() << "\""
                  << ","    << (*bin)->Flux()
                  << ","    << (*bin)->FluxPErr()
                  << ","    << (*bin)->FluxNErr()
                  << ",\""  << spectrum_iter->second->FluxUnits()   << "\""
                  << ");";

              SLOG("NuVld", pINFO)
                      << "SQL string to be sent to the DBase: "
                                                     << sql_string.str().c_str();

              _sql_server->Query( sql_string.str().c_str() );

       } // loop ober spectrum bins
  } // loop over spectra

  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadReferences(
       const XmlExperimentInfo & exp_info, const XmlMeasurementHeader & m_header) const
{
  const vector<XmlCitation *> & refs = m_header.GetRefs();

  int iref=0;
  vector<XmlCitation *>::const_iterator ref_iter;
  for(ref_iter = refs.begin(); ref_iter != refs.end(); ++ref_iter) {

     ostringstream sql_string;

     sql_string << "INSERT INTO REFERENCE VALUES "
                << "(\""  << exp_info.Name()         << "\""
                << ",\""  << m_header.Tag()          << "\""
                << ","    << iref++
                << ",\""  << (*ref_iter)->Author()   << "\""
                << ",\""  << (*ref_iter)->Journal()  << "\""
                << ","    << (*ref_iter)->Year()
                << ");";

     SLOG("NuVld", pINFO)
          << "SQL string to be sent to the DBase: " << sql_string.str().c_str();

     _sql_server->Query( sql_string.str().c_str() );

  } //loop over references

  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadMeasHdr(
 const XmlExperimentInfo & exp_info, const XmlMeasurementHeader & m_header) const
{
  ostringstream sql_string;

  sql_string << "INSERT INTO MEASUREMENT_HEADER VALUES "
                << "(\""  << exp_info.Name()           << "\""
                << ",\""  << m_header.Tag()            << "\""
                << ",\""  << m_header.Observable()     << "\""
                << ",\""  << m_header.Target()         << "\""
                << ",\""  << m_header.Reaction()       << "\""
                << ",\""  << m_header.A()              << "\""
                << ",\""  << m_header.Exposure()       << "\""
                << ",\""  << m_header.ExposureUnits()  << "\""
                << ",\""  << m_header.DataSource()     << "\""
                << ","    << m_header.NPoints()
                << ","    << m_header.ErrorStatus()
                << ","    << m_header.GetRefs().size()
                << ",\""  << m_header.Comment()        << "\""
                << ");";

  SLOG("NuVld", pINFO)
     << "\nSQL string to be sent to the DBase: " << sql_string.str().c_str();

  _sql_server->Query( sql_string.str().c_str() );

  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadNuXSec(
    const XmlExperimentInfo & exp_info, const XmlMeasurement * meas) const
{
  const vector<XmlRecordBase *> & d_points = meas->GetDataPoints();
  const XmlMeasurementHeader &    m_header = meas->GetXmlMeasurementHeader();

  vector<XmlRecordBase *>::const_iterator point_iter = d_points.begin();
  for( ; point_iter != d_points.end(); ++point_iter){

     ostringstream sql_string;

     sql_string << "INSERT INTO CROSS_SECTION VALUES "
                << "(\""  << exp_info.Name()               << "\""
                << ",\""  << m_header.Tag()                << "\""
                << ","    << (*point_iter)->Get("xsec")
                << ","    << (*point_iter)->Get("stat_err+")
                << ","    << (*point_iter)->Get("stat_err-")
                << ","    << (*point_iter)->Get("syst_err+")
                << ","    << (*point_iter)->Get("syst_err-")
                << ",\""  << (*point_iter)->Get("xsec_units")     << "\""
                << ",\""  << (*point_iter)->Get("xsec_norm")      << "\""
                << ",\""  << (*point_iter)->Get("stat_err_value") << "\""
                << ",\""  << (*point_iter)->Get("syst_err_value") << "\""
                << ","    << (*point_iter)->Get("E")
                << ","    << (*point_iter)->Get("E_min")
                << ","    << (*point_iter)->Get("E_max")
                << ",\""  << (*point_iter)->Get("E_units")  << "\""
                << ",\""  << (*point_iter)->Get("E_frame")  << "\""
                << ");";

     SLOG("NuVld", pINFO)
          << "SQL string to be sent to the DBase: " << sql_string.str().c_str();

     _sql_server->Query( sql_string.str().c_str() );
  }
  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadElDiffXSec(
    const XmlExperimentInfo & exp_info, const XmlMeasurement * meas) const
{
  const vector<XmlRecordBase *> & d_points = meas->GetDataPoints();
  const XmlMeasurementHeader &    m_header = meas->GetXmlMeasurementHeader();

  vector<XmlRecordBase *>::const_iterator point_iter = d_points.begin();

  for( ; point_iter != d_points.end(); ++point_iter){

     ostringstream sql_string;

     sql_string << "INSERT INTO E_DIFF_CROSS_SECTION VALUES "
                << "(\""  << exp_info.Name()               << "\""
                << ",\""  << m_header.Tag()                << "\""
                << ","    << (*point_iter)->Get("Sigma")
                << ",\""  << (*point_iter)->Get("Sigma_units") << "\""
                << ","    << (*point_iter)->Get("dSigma")
                << ","    << (*point_iter)->Get("E")
                << ",\""  << (*point_iter)->Get("E_units")     << "\""
                << ","    << (*point_iter)->Get("Ep")
                << ",\""  << (*point_iter)->Get("Ep_units")    << "\""
                << ","    << (*point_iter)->Get("Theta")
                << ",\""  << (*point_iter)->Get("Theta_units") << "\""
                << ","    << (*point_iter)->Get("Q2")
                << ",\""  << (*point_iter)->Get("Q2_units")    << "\""
                << ","    << (*point_iter)->Get("W2")
                << ",\""  << (*point_iter)->Get("W2_units")    << "\""
                << ","    << (*point_iter)->Get("Nu")
                << ",\""  << (*point_iter)->Get("Nu_units")    << "\""
                << ","    << (*point_iter)->Get("Epsilon")
                << ","    << (*point_iter)->Get("Gamma")
                << ","    << (*point_iter)->Get("x")
                << ");";

     SLOG("NuVld", pINFO)
          << "SQL string to be sent to the DBase: " << sql_string.str().c_str();

     _sql_server->Query( sql_string.str().c_str() );
  }
  return eDbu_OK;
}
//_______________________________________________________________________________
DBStatus_t DBXmlUploader::UploadSF(
   const XmlExperimentInfo & exp_info, const XmlMeasurement * meas) const
{
  const vector<XmlRecordBase *> & d_points = meas->GetDataPoints();
  const XmlMeasurementHeader &    m_header = meas->GetXmlMeasurementHeader();

  vector<XmlRecordBase *>::const_iterator point_iter = d_points.begin();

  for( ; point_iter != d_points.end(); ++point_iter){

     ostringstream sql_string;

     sql_string << "INSERT INTO STRUCTURE_FUNCTION VALUES "
                << "(\""  << exp_info.Name()            << "\""
                << ",\""  << m_header.Tag()             << "\""
                << ","    << (*point_iter)->Get("sf")
                << ",\""  << (*point_iter)->Get("sf_p") << "\""
                << ",\""  << (*point_iter)->Get("sf_R") << "\""
                << ","    << (*point_iter)->Get("x")
                << ","    << (*point_iter)->Get("Q2")
                << ","    << (*point_iter)->Get("stat_err+")
                << ","    << (*point_iter)->Get("stat_err-")
                << ","    << (*point_iter)->Get("syst_err+")
                << ","    << (*point_iter)->Get("syst_err-")
                << ");";

     SLOG("NuVld", pINFO)
          << "SQL string to be sent to the DBase: " << sql_string.str().c_str();

     _sql_server->Query( sql_string.str().c_str() );
  }
  return eDbu_OK;
}
//_______________________________________________________________________________

