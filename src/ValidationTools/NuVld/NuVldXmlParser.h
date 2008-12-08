//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldXmlParser

\brief    A libxml2 based parser for the NuValidator XML data files

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _NUVLD_XML_PARSER_H_
#define _NUVLD_XML_PARSER_H_

#include <string>
#include <map>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include "ValidationTools/NuVld/XmlObservable.h" 
#include "ValidationTools/NuVld/XmlDataSet.h"
#include "ValidationTools/NuVld/XmlExperimentMeasurements.h"
#include "ValidationTools/NuVld/XmlMeasurement.h"
#include "ValidationTools/NuVld/XmlExperimentInfo.h"
#include "ValidationTools/NuVld/XmlCitation.h"
#include "ValidationTools/NuVld/XmlRecord.h"
#include "ValidationTools/NuVld/ParserStatus.h"

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
  XmlExperimentMeasurements * ParseExperiment        (xmlNodePtr xml_cur, string name);
  XmlExperimentInfo *         ParseXmlExperimentInfo    (xmlNodePtr xml_cur);
  XmlBeamFluxSpectrum *       ParseXmlBeamFluxSpectrum  (xmlNodePtr xml_cur);
  XmlBeamFluxBin *            ParseXmlBeamFluxBin       (xmlNodePtr xml_cur);
  XmlMeasurement *            ParseXmlMeasurement       (xmlNodePtr xml_cur, XmlObservable_t obs);
  XmlMeasurementHeader *      ParseXmlMeasurementHeader (xmlNodePtr xml_cur);
  XmlCitation *               ParseReference         (xmlNodePtr xml_cur);
  XmlRecordBase *             ParsePoint             (xmlNodePtr xml_cur, XmlObservable_t obs);
  XmlRecordBase *             CreateNewXmlRecord        (XmlObservable_t obs);
      
  xmlDocPtr                _xml_doc;
  string                   _xml_filename;   
  XmlDataSet *             _data;
  int                      _measurement_tag;
  XmlParserStatus_t        _parser_status;
};

} // nuvld namespace
} // genie namespace

#endif // _PARSER_H_
