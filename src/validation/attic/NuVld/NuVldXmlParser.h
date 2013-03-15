//____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldXmlParser

\brief    A libxml2-based parser for the NuValidator XML data files

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  August 2003 

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUVLD_XML_PARSER_H_
#define _NUVLD_XML_PARSER_H_

#include <string>
#include <map>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include "Conventions/XmlParserStatus.h"
#include "ValidationTools/NuVld/XmlObservable.h" 
#include "ValidationTools/NuVld/XmlDataSet.h"
#include "ValidationTools/NuVld/XmlExperimentMeasurements.h"
#include "ValidationTools/NuVld/XmlMeasurement.h"
#include "ValidationTools/NuVld/XmlExperimentInfo.h"
#include "ValidationTools/NuVld/XmlCitation.h"
#include "ValidationTools/NuVld/XmlRecord.h"

using std::string;
using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace genie {
namespace nuvld {
  
class NuVldXmlParser {

public:
  NuVldXmlParser();
 ~NuVldXmlParser();

  //! parse the input XML document
  void ParseXmlDocument(const char * filename);

  //! return the parser status  
  XmlParserStatus_t GetXmlParsingStatus(void) const { return fXmlPStatus; }
  
  //! get the loaded data set
  const XmlDataSet & GetDataSet(void) const;

private:

  XmlParserStatus_t           Verify                 (void);
  void                        FillDataSet            (void);
  XmlExperimentMeasurements * ParseExperiment        (xmlNodePtr xml_cur, string name);
  XmlExperimentInfo *         ParseExperimentInfo    (xmlNodePtr xml_cur);
  XmlBeamFluxSpectrum *       ParseBeamFluxSpectrum  (xmlNodePtr xml_cur);
  XmlBeamFluxBin *            ParseBeamFluxBin       (xmlNodePtr xml_cur);
  XmlMeasurement *            ParseMeasurement       (xmlNodePtr xml_cur, XmlObservable_t obs, int id);
  XmlMeasurementHeader *      ParseMeasurementHeader (xmlNodePtr xml_cur);
  XmlCitation *               ParseReference         (xmlNodePtr xml_cur);
  XmlRecordBase *             ParsePoint             (xmlNodePtr xml_cur, XmlObservable_t obs);
  XmlRecordBase *             CreateNewXmlRecord     (XmlObservable_t obs);
      
  xmlDocPtr         fXmlDoc;       ///< parsed XML doc
  string            fXmlFilename;  ///< input XML file name
  XmlDataSet *      fDataSet;      ///< parsed data set
  XmlParserStatus_t fXmlPStatus;   ///< parser status
};

} // nuvld namespace
} // genie namespace

#endif // _PARSER_H_
