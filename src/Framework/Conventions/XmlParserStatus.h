//____________________________________________________________________________
/*!

\class    genie::XmlParserStatus

\brief    Encapsulates an XML document parsing status

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 4, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//_______________________________________________________________________________

#ifndef _XML_PARSER_STATUS_H_
#define _XML_PARSER_STATUS_H_

namespace genie {

typedef enum EXmlParseStatus {

  kXmlUndefined   = -1,
  kXmlOK          =  0,
  kXmlNotParsed   =  1,
  kXmlEmpty       =  2,
  kXmlInvalidRoot =  3         

} XmlParserStatus_t;


class XmlParserStatus {

  public:
     static const char * AsString(XmlParserStatus_t status) {
       switch(status) {
         case kXmlUndefined:   return "Undefined state";                       break;
         case kXmlOK:          return "XML document succesfully parsed";       break;
         case kXmlNotParsed:   return "XML document parsing failed";           break;
         case kXmlEmpty:       return "XML document is empty";                 break;
         case kXmlInvalidRoot: return "XML document has invalid root element"; break;
         default:              return "unrecognized XMLParseStatus_t enum";    break;
       }
       return "unrecognized XMLParseStatus_t enum";  
     }
};

}       // namespace
#endif  // _XML_PARSER_STATUS_H_

