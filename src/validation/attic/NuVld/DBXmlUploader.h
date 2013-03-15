//____________________________________________________________________________
/*!

\class    genie::nuvld::DBXmlUploader

\brief    Utility class used by the DBI to upload XML data to the RDBMS

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Jan, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DBASE_UPLOADER_H_
#define _DBASE_UPLOADER_H_

#include <TSQLServer.h>

#include "ValidationTools/NuVld/XmlDataSet.h"
#include "ValidationTools/NuVld/XmlExperimentInfo.h"
#include "ValidationTools/NuVld/XmlMeasurement.h"
#include "ValidationTools/NuVld/XmlExperimentMeasurements.h"
#include "ValidationTools/NuVld/DBStatus.h"

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
  
  DBStatus_t UploadExpInfo     (const XmlExperimentInfo & expi) const; 
  DBStatus_t UploadFluxSpectra (const XmlExperimentInfo & expi,
                                const map<string, XmlBeamFluxSpectrum *> & spectra) const;
  DBStatus_t UploadReferences  (const XmlExperimentInfo & expi,
                                const XmlMeasurementHeader & mhdr) const;
  DBStatus_t UploadMeasHdr     (const XmlExperimentInfo & expi,
                                const XmlMeasurementHeader & mhdr) const;
  DBStatus_t UploadNuXSec      (const XmlExperimentInfo & expi, const XmlMeasurement * meas) const;
  DBStatus_t UploadElDiffXSec  (const XmlExperimentInfo & expi, const XmlMeasurement * meas) const;
  DBStatus_t UploadSF          (const XmlExperimentInfo & expi, const XmlMeasurement * meas) const;

  TSQLServer * _sql_server;
};

} // nuvld namespace
} // genie namespace

#endif
