//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBXmlUploader

\brief    Utility class used by the DBI to upload XML data to the RDBMS

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _DBASE_UPLOADER_H_
#define _DBASE_UPLOADER_H_

#include <TSQLServer.h>

#include "XmlParser/XmlDataSet.h"
#include "XmlParser/ExperimentInfo.h"
#include "XmlParser/Measurement.h"
#include "XmlParser/ExperimentMeasurements.h"
#include "DBUtils/DBStatus.h"

using namespace std;
using namespace genie::nuvld;

namespace genie {
namespace nuvld {

class DBXmlUploader 
{
  friend class DBI;

private:

  DBXmlUploader(TSQLServer * server);
  virtual ~DBXmlUploader();

  DBStatus_t Upload            (const XmlDataSet & data) const;
  
  DBStatus_t UploadExpInfo     (const ExperimentInfo & expi) const; 
  DBStatus_t UploadFluxSpectra (const ExperimentInfo & expi,
                                const map<string, BeamFluxSpectrum *> & spectra) const;
  DBStatus_t UploadReferences  (const ExperimentInfo & expi,
                                const MeasurementHeader & mhdr) const;
  DBStatus_t UploadMeasHdr     (const ExperimentInfo & expi,
                                const MeasurementHeader & mhdr) const;
  DBStatus_t UploadNuXSec      (const ExperimentInfo & expi, const Measurement * meas) const;
  DBStatus_t UploadElDiffXSec  (const ExperimentInfo & expi, const Measurement * meas) const;
  DBStatus_t UploadSF          (const ExperimentInfo & expi, const Measurement * meas) const;

  TSQLServer *       _sql_server;
};

} // nuvld namespace
} // genie namespace

#endif
