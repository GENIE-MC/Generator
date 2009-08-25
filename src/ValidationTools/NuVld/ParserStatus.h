//_____________________________________________________________________________
/*!

\class    genie::nuvld::ParserStatus

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#ifndef _PARSER_STATUS_H_
#define _PARSER_STATUS_H_

namespace genie {
namespace nuvld {

typedef enum EXmlParseStatus {

   eXml_OK           = 0,
   eXml_NOT_PARSED   = 1,
   eXml_EMPTY        = 2,
   eXml_INVALID_ROOT = 3         

} XmlParserStatus_t;

class ParserStatus {

  public:

     static const char * AsString(XmlParserStatus_t status) {
       switch(status) {
         case eXml_OK:           return "XML document succesfully parsed";       break;
         case eXml_NOT_PARSED:   return "XML document parsing failed";           break;
         case eXml_EMPTY:        return "XML document is empty";                 break;
         case eXml_INVALID_ROOT: return "XML document has invalid root element"; break;
         default:                return "unrecognized XML parser status";        break;
       }
       return "unrecognized XML parser status";  
     }
};

} // nuvld namespace
} // genie namespace

#endif // _PARSER_STATUS_H_

