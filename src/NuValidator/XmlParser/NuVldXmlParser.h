//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldXmlParser

\brief    A libxml2 based parser for the NuValidator XML data files

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _NUVLD_XML_PARSER_H_
#define _NUVLD_XML_PARSER_H_

#include <string>
#include <map>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include "XmlParser/Observable.h" 
#include "XmlParser/XmlDataSet.h"
#include "XmlParser/ExperimentMeasurements.h"
#include "XmlParser/Measurement.h"
#include "XmlParser/ExperimentInfo.h"
#include "XmlParser/Citation.h"
#include "XmlParser/Record.h"
#include "XmlParser/ParserStatus.h"

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

  void ParseXmlDocument(const char * filename);
  
  XmlParserStatus_t GetXmlParsingStatus(void) const { return _parser_status; }
  
  const XmlDataSet & GetDataSet(void) const;

private:

  XmlParserStatus_t VerifyParsing(void);

  void                     FillDataSet            (void);
  ExperimentMeasurements * ParseExperiment        (xmlNodePtr xml_cur, string name);
  ExperimentInfo *         ParseExperimentInfo    (xmlNodePtr xml_cur);
  BeamFluxSpectrum *       ParseBeamFluxSpectrum  (xmlNodePtr xml_cur);
  BeamFluxBin *            ParseBeamFluxBin       (xmlNodePtr xml_cur);
  Measurement *            ParseMeasurement       (xmlNodePtr xml_cur, Observable_t obs);
  MeasurementHeader *      ParseMeasurementHeader (xmlNodePtr xml_cur);
  Citation *               ParseReference         (xmlNodePtr xml_cur);
  RecordBase *             ParsePoint             (xmlNodePtr xml_cur, Observable_t obs);
  RecordBase *             CreateNewRecord        (Observable_t obs);
      
  xmlDocPtr                _xml_doc;
  string                   _xml_filename;   
  XmlDataSet *             _data;
  int                      _measurement_tag;
  XmlParserStatus_t        _parser_status;
};

} // nuvld namespace
} // genie namespace

#endif // _PARSER_H_
